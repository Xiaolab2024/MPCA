library(tidyverse)
library(glmnet)
library(ggpubr)
library(data.table)
library(matrixStats)
library(matrixTests)
library(ggrepel)
library(Hmisc)
library(fdrtool)

### read in data frame
df<-read.csv("(0)matrix_protein_metabolite.csv",head=T,stringsAsFactors = F)

## transpose
df1<-df%>%unite("ID",c("analyte","type"),sep="^")
dft<-as.data.frame(t(df1[,-1]))
colnames(dft)<-df1$ID

matrix<-as.matrix(dft)
matrix[is.nan(matrix)]<-NA
matrix[is.infinite(matrix)] <- NA

##pearson correlation
p2p_pearson<-rcorr(matrix, type="pearson")

pmat<-p2p_pearson$P

cormat<-p2p_pearson$r

nmat<-p2p_pearson$n

flattenCorrMatrix <- function(cormat, pmat,nmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut],
    n = nmat[ut],
    stringsAsFactors = FALSE
  )
}


#flatten the pearson table and remove correlations calculated with <50 paired measurements
flat_pearson<-flattenCorrMatrix(cormat, pmat,nmat)
flat_pearson_filtered<-flat_pearson%>%
  filter(n>50)

flat_pearson_filtered1<-flat_pearson_filtered%>%
  separate(row,c("row","analyte1"),sep="\\^")%>%
  separate(column,c("column","analyte2"),sep="\\^")

flat_pearson_filtered2<-flat_pearson_filtered1%>%
  filter(analyte1 != analyte2)

#BH p value adjust
flat_pearson_filtered2$"p_adj"<-p.adjust(flat_pearson_filtered2$p, method = "BH")
flat_pearson_filtered3<-flat_pearson_filtered2%>%
  filter(p_adj< 0.05)
