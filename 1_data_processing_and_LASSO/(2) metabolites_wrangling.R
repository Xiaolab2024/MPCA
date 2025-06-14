library(tidyverse)
library(matrixStats)
library(matrixTests)
library(ggpubr)


#metabolites measured under negative mode
neg<-read.csv("BAT_neg.csv",header=T,stringsAsFactors = F)%>%arrange(Sample)

#remove all 0 columns
neg <- neg[, base::colSums(neg != 0) > 0]
neg<-neg[- grep("STD", neg$Sample),]
neg<-neg[- grep("poolid", neg$Sample),]

# remove metabolites that are 10000 fold lower than avg. internal standard
avg_neg<-mean(base::colSums(as.data.frame(neg$D4.thymine))/nrow(neg)+
                base::colSums(as.data.frame(neg$Glycocholate.d4))/nrow(neg)+
                base::colSums(as.data.frame(neg$Inosine.15N4))/nrow(neg))

neg_mat <- neg%>%select(-Sample)%>%select_if(~median(., na.rm = TRUE) >= (avg_neg/10000))
Sample<-neg$Sample
Sample<-as.data.frame(Sample)

neg_filtered<-cbind(Sample,neg_mat)%>%dplyr::rename("D4.thymine_neg"="D4.thymine",
                                                    "Glycocholate.d4_neg"="Glycocholate.d4",
                                                    "Inosine.15N4_neg"="Inosine.15N4")

# get the sample-to-bridge(pool) ratios
neg_filtered<-neg_filtered%>%filter(Sample!=	"X13618_pool_neg" )
negt<-as.data.frame(t(neg_filtered[,2:ncol(neg_filtered)]))
rownames(negt)<-colnames(neg_filtered[,2:ncol(neg_filtered)])
colnames(negt)<-neg_filtered$Sample

negt_pools<-negt%>%select(ends_with("_pool_neg"))
negt_mocks<-negt%>%select(ends_with("_poolmock_neg"))
negt_ratio<-negt%>%select(-ends_with("_pool_neg"))



negt_ratio[,1:11]<-negt[,1:11]/negt[,12]## has mock
negt_ratio[,12:21]<-negt[,13:22]/negt[,23]
negt_ratio[,22:31]<-negt[,24:33]/negt[,34]
negt_ratio[,32:41]<-negt[,35:44]/negt[,45]
negt_ratio[,42:51]<-negt[,46:55]/negt[,56]
negt_ratio[,52:62]<-negt[,57:67]/negt[,68]## has mock
negt_ratio[,63:72]<-negt[,69:78]/negt[,79]
negt_ratio[,73:82]<-negt[,80:89]/negt[,90]
negt_ratio[,83:92]<-negt[,91:100]/negt[,101]
negt_ratio[,93:102]<-negt[,102:111]/negt[,112]
negt_ratio[,103:113]<-negt[,113:123]/negt[,124]## has mock
negt_ratio[,114:123]<-negt[,125:134]/negt[,135]
negt_ratio[,124:133]<-negt[,136:145]/negt[,146]
negt_ratio[,134:143]<-negt[,147:156]/negt[,157]
negt_ratio[,144:153]<-negt[,158:167]/negt[,168]
negt_ratio[,154:164]<-negt[,169:179]/negt[,180]## has mock
negt_ratio[,165:171]<-negt[,181:187]/negt[,188]## has fewer samples

neg_ratio<-as.data.frame(t(negt_ratio))
neg_ratio$Sample<-rownames(neg_ratio)
neg_ratio<-neg_ratio%>%
  mutate(Mouse="M")%>%
  separate("Sample",c("a","ID","b"),sep="_",remove=F)%>%
  select(-a,-b)
neg_ratio_null<-neg_ratio
rownames(neg_ratio_null)<-NULL
neg_ratio$number<-rownames(neg_ratio_null)

neg_ratio<-neg_ratio%>%
  mutate(ID=case_when(
    ID== "poolmock" ~ paste0(ID,number),
    TRUE ~ ID
  ))

neg_ratio<-neg_ratio%>%
  unite("Sample",c("Mouse","ID"),sep="_",remove=T)%>%
  select(-number)
######## write_tsv(neg_ratio,"bat_neg_ratio.tsv")


######## calibrate by internal standard ######
neg_ratio1<-neg_ratio%>%mutate(norm_facotr=(D4.thymine_neg+Glycocholate.d4_neg+Inosine.15N4_neg)/3)
neg_ratio2<-cbind(neg_ratio1[,1:(ncol(neg_ratio1)-2)]/neg_ratio1[,ncol(neg_ratio1)],as.data.frame(neg_ratio1$Sample)%>%dplyr::rename('Sample'='neg_ratio1$Sample'))
######## write_tsv(neg_ratio2,"bat_neg_ratio_injcali.tsv")




##########metabolites measured under positive mode######
pos<-read.csv("BAT_pos.csv",header=T,stringsAsFactors = F)%>%arrange(Sample)

#remove all 0 columns
pos <- pos[, base::colSums(pos != 0) > 0]
#pos<-pos[- grep("mock", pos$Sample),]
pos<-pos[- grep("STD", pos$Sample),]
pos<-pos[- grep("poolid", pos$Sample),]

# remove metabolites that are 10000 fold lower than avg. internal standard
avg_pos<-mean(  base::colSums(as.data.frame(pos$Glycocholate.d4))/nrow(pos)+
                  base::colSums(as.data.frame(pos$Inosine.15N4))/nrow(pos))

pos_mat <- pos%>%select(-Sample)%>%select_if(~median(., na.rm = TRUE) >= (avg_pos/100))
Sample<-pos$Sample
Sample<-as.data.frame(Sample)

pos_filtered<-cbind(Sample,pos_mat)%>%dplyr::rename("Glycocholate.d4_pos"="Glycocholate.d4",
                                                    "Inosine.15N4_pos"="Inosine.15N4")

# get the sample-to-bridge(pool) ratios
pos_filtered<-pos_filtered%>%filter(Sample!=	"X13408_pool_pos" )
post<-as.data.frame(t(pos_filtered[,2:ncol(pos_filtered)]))
rownames(post)<-colnames(pos_filtered[,2:ncol(pos_filtered)])
colnames(post)<-pos_filtered$Sample

post_pools<-post%>%select(ends_with("_pool_pos"))
post_mocks<-post%>%select(ends_with("_poolmock_pos"))
post_ratio<-post%>%select(-ends_with("_pool_pos"))


post_ratio[,1:11]<-post[,1:11]/post[,12]## has mock
post_ratio[,12:21]<-post[,13:22]/post[,23]
post_ratio[,22:31]<-post[,24:33]/post[,34]
post_ratio[,32:41]<-post[,35:44]/post[,45]
post_ratio[,42:51]<-post[,46:55]/post[,56]
post_ratio[,52:62]<-post[,57:67]/post[,68]## has mock
post_ratio[,63:72]<-post[,69:78]/post[,79]
post_ratio[,73:82]<-post[,80:89]/post[,90]
post_ratio[,83:92]<-post[,91:100]/post[,101]
post_ratio[,93:102]<-post[,102:111]/post[,112]
post_ratio[,103:113]<-post[,113:123]/post[,124]## has mock
post_ratio[,114:123]<-post[,125:134]/post[,135]
post_ratio[,124:133]<-post[,136:145]/post[,146]
post_ratio[,134:143]<-post[,147:156]/post[,157]
post_ratio[,144:153]<-post[,158:167]/post[,168]
post_ratio[,154:164]<-post[,169:179]/post[,180]## has mock
post_ratio[,165:171]<-post[,181:187]/post[,188]## has fewer samples

pos_ratio<-as.data.frame(t(post_ratio))
pos_ratio$Sample<-rownames(pos_ratio)
pos_ratio<-pos_ratio%>%
  mutate(Mouse="M")%>%
  separate("Sample",c("a","ID","b"),sep="_",remove=F)%>%
  select(-a,-b)

pos_ratio_null<-pos_ratio
rownames(pos_ratio_null)<-NULL
pos_ratio$number<-rownames(pos_ratio_null)

pos_ratio<-pos_ratio%>%
  mutate(ID=case_when(
    ID== "poolmock" ~ paste0(ID,number),
    TRUE ~ ID
  ))

pos_ratio<-pos_ratio%>%
  unite("Sample",c("Mouse","ID"),sep="_",remove=T)%>%
  select(-number)
######## write_tsv(pos_ratio,"bat_pos_ratio.tsv")

######## calibrate by internal standard ######
pos_ratio1<-pos_ratio%>%mutate(norm_facotr=(Glycocholate.d4_pos+Inosine.15N4_pos)/2)
pos_ratio2<-cbind(pos_ratio1[,1:(ncol(pos_ratio1)-2)]/pos_ratio1[,ncol(pos_ratio1)],as.data.frame(pos_ratio1$Sample)%>%dplyr::rename('Sample'='pos_ratio1$Sample'))
######## write_tsv(pos_ratio2,"bat_pos_ratio_injcali.tsv")
