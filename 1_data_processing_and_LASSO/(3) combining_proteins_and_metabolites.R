library(tidyverse)

#read in proteins
dfp<-read.delim("protein.tsv",header = T,stringsAsFactors = F)
dfp<-dfp[,-c(3:18)]
dfp1<-dfp%>%
  unite(analyte,c("Protein.Id","Gene.Symbol"),sep="$")%>%
  mutate("type"="protein")


#read in metabolites
dfm<-read.csv("(0) metabolite.csv",header = T,stringsAsFactors = F)

#combine data
df_all<- bind_rows(dfp1, dfm)
df_all<-df_all[,c(1,ncol(df_all),2:(ncol(df_all)-1))]

df_all<-as.matrix(df_all)

#remove NA and inf
df_all[is.nan(df_all)]<-NA
df_all[is.infinite(df_all)] <- NA

#final dataset
df_all<-as.data.frame(df_all)


