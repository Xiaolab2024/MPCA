library(tidyverse) 
library(glmnet)
library(forcats)
library(ggpubr)
library(ggridges)

#read in data
bat_all<- read.delim("matrix_protein_metabolite.tsv",header=TRUE, stringsAsFactors = F)

#get proteins
bat_protein<-bat_all%>%
  filter(type=="protein")%>%
  select(-type)
#get metabolites
depvar0<- bat_all%>%
  filter(type=="metabolite")%>%
  select(-type)
# format the table
depvar<-as.data.frame(t(depvar0[,2:ncol(depvar0)]))
colnames(depvar)<-depvar0$analyte
rownames(depvar)<-colnames(depvar0[,2:ncol(depvar0)])
depvar$Sample_ID1<-rownames(depvar)

#replace possible inf with NA
bat_protein<-do.call(data.frame,lapply(bat_protein, function(x) replace(x, is.infinite(x),NA)))

# remove entries with percent NA values in more than 50 samples
bat_protein<-bat_protein[rowMeans(is.na(bat_protein[,2:ncol(bat_protein)]))<.306,]
col_names<-bat_protein$analyte
row_names<-colnames(bat_protein[,2:ncol(bat_protein)])
df<-as.data.frame(t(bat_protein[,2:ncol(bat_protein)]))
colnames(df)<-col_names
rownames(df)<-row_names
df$Sample_ID1<-rownames(df)
indep_var_scaled<-df ## no scaling is used despite the variable name
protein_scaled <- indep_var_scaled
db_select<-depvar

#make sure tables are sorted before run next line!!!
protein_scaled<-left_join(protein_scaled,db_select,by="Sample_ID1")%>%select(-Sample_ID1)

####lasso blow###
indep_var_scaled_noNA<-indep_var_scaled[ , colSums(is.na(indep_var_scaled)) == 0]%>%select(-Sample_ID1)
indep_var_scaled_1<-indep_var_scaled%>%select(-Sample_ID1)
indep_var_scaled_noNA_scale<- data.frame(lapply(indep_var_scaled_noNA, function(x) scale(x, center = FALSE, scale = max(x, na.rm = TRUE)/100)))

#for loops for lasso and plots##

for (i in colnames(depvar)){

list<-as.matrix(depvar%>%dplyr::select(i))

#non-scaled

if(any(is.na(list))==T){
  next
} else{

set.seed=100000000
lasso_fat_mass <-cv.glmnet(as.matrix(indep_var_scaled_noNA), list, weights = NULL, offset = NULL, lambda = NULL,
                           type.measure = c("default", "mse", "deviance", "class", "auc", "mae","C"), 
                           nfolds = 10, foldid = NULL, alignment = c("lambda",  "fraction"), 
                           grouped = TRUE, keep = FALSE, parallel = FALSE,
                           gamma = c(0, 0.25, 0.5, 0.75, 1), relax = FALSE, trace.it = 1)


plot(lasso_fat_mass)
coef(lasso_fat_mass)
lasso_fat_mass

coef_lasso<-as.matrix(coef(lasso_fat_mass)[c(2:nrow(coef(lasso_fat_mass))),])
coef_matrix<-as.data.frame(coef_lasso)
coef_matrix$metabolites<-rownames(coef_matrix)
coef_variable<-coef_matrix[coef_matrix$V1!=0,]

indep_var_scaled_noNA_select<-indep_var_scaled_noNA%>%select(coef_variable$metabolites)
indep_var_scaled_noNA_select$fat_mass<-depvar$fat_mass


##plot lasso predictors##
if(nrow(coef_variable)==0){
  next
} else {
hits_plot<-coef_variable %>%
  mutate(name = fct_reorder(metabolites, V1)) %>%
  ggplot( aes(x=name, y=V1)) +
  ylab(paste0(i,"_coefficient"))+
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
  coord_flip() +
  xlab("") +
  theme_bw()


##pearson r of prediction correlating predicted value to true value##
protein_scaled$lasso_predict<-predict(lasso_fat_mass, as.matrix(indep_var_scaled_noNA))
prediction_plot<-ggscatter(protein_scaled, x = "lasso_predict", y = i, 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "lasso_prediction", ylab = paste0(i,"_abundance"))

coef_matrix1<-coef_matrix%>%
  dplyr::select(metabolites, V1)%>%
  dplyr::rename(!!i:=V1)

i1<-gsub("\\/", "_or_", i)
ggsave(paste0("meta_",i1,"_hitsplot.pdf"), plot=hits_plot,width = 10, height = 10, units = "in")
ggsave(paste0("meta_",i1,"_predictionplot.pdf"),plot=prediction_plot, width = 5, height = 5, units = "in")

if(i==colnames(depvar)[[1]]){
  
  coef_final<-coef_matrix1
  
  
}else{
  coef_final<-left_join(coef_final,coef_matrix1,by="metabolites")
}}}

print(i)
}

