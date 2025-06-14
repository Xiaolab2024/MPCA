library(tidyverse)

##read in data files
names<-read.csv("(1) channel_mouse_matches.csv",header=T,stringsAsFactors = F)%>%
  mutate(a="M")%>%
  unite("Sample_ID","a":"Sample_ID",remove=TRUE)
df<-read.delim("(0) protein_quant.tsv", header=T,stringsAsFactors = F)
df<-df[-grep("##|contaminant",df$Protein.Id),]

##divide_by_bridge_and_log2,change c(1:x) based on number of sets
for(i in c(1:12)){
  df_samp<-log2(
           df[(19+(i-1)*16):(33+(i-1)*16)]/df[,(34+(i-1)*16)]
             )
        
  if(i==1){
    final<-cbind(df[1:18],df_samp)
  }
    else {
      final<-cbind(final,df_samp)
    }
}

#match TMT channels to mouse ID
final<-final%>%
  select(-L.rq_126_sn.sum	,-L.rq_127n_sn.sum,-L.rq_127c_sn.sum,-L.rq_128n_sn.sum,-L.rq_128c_sn.sum,
         -L.rq_129n_sn.sum,-L.rq_129c_sn.sum,-L.rq_130c_sn.sum,-L.rq_131_sn.sum,-L.rq_131c_sn.sum,
         -L.rq_132n_sn.sum,-L.rq_133n_sn.sum,-L.rq_133c_sn.sum)
tfinal<-as.data.frame(t(final[,19:ncol(final)]))
tfinal$channel<-rownames(tfinal)
tfinal<-tfinal[,c(ncol(tfinal),1:(ncol(tfinal)-1))]
join<-left_join(tfinal,names,by="channel")
join<-join[,c(ncol(join),1:(ncol(join)-1))]
row.names(join)<-join$Sample_ID
join<-join[,3:ncol(join)]
final2<-as.data.frame(t(join))

#final protein table
final3<-cbind(final[,1:18],final2)
