# This script is to count the conservation of GLK target genes in five species.
# The Orthogroups.tsv result is got by OrthoFinder.
data<-read.table('../../data/Species_specific_analysis_data/Orthogroups.tsv',sep=':')
# Calculate the conservation and add orthogroup id on the gene table.
conserve<-NULL
for (i in 1:nrow(data)){
conserve[i]<-sum(grepl('AT',data[i,2]),grepl('Niben',data[i,2]),grepl('LOC',data[i,2]),grepl('Soly',data[i,2]),grepl('Zm',data[i,2]))
}
data<-data.frame(data,conserve=conserve)
data<-as.matrix(data)
gene_group<-apply(data,1,function(x) {
gene<-unlist(strsplit(x[2],' '))
gene<-gene[gene!='']
result<-data.frame(Ortho_id=rep(x[1],length(gene)),gene=gene,conserve=rep(x[3],length(gene)))
return(result)
}
)
conserve_count<-Reduce(rbind,gene_group)
# Write the table by the species.
AT<-conserve_count[grep('AT',conserve_count[,2]),]
NB<-conserve_count[grep('Niben',conserve_count[,2]),]
SL<-conserve_count[grep('Soly',conserve_count[,2]),]
OS<-conserve_count[grep('LOC',conserve_count[,2]),]
ZM<-conserve_count[grep('Zm',conserve_count[,2]),]

# Combine the ChIP-seq information.
# All these tables are in Supplementary S3.
AT_table<-read.table('~/glk/At_960.csv',sep='\t',header=T,quote='',fill=T)
NB_table<-read.table('~/glk/Nben_956.csv',sep='\t',header=T,quote='',fill=T)
SL_table<-read.table('~/glk/tomato_1286.csv',sep='\t',header=T,quote='',fill=T)
OS_table<-read.table('~/glk/rice_332.csv',sep='\t',header=T,quote='',fill=T)
ZM_table<-read.table('~/glk/maize_1089.csv',sep='\t',header=T,quote='',fill=T)
AT_table<-merge(AT_table,AT,by='gene',all.x=T)
NB_table<-merge(NB_table,NB,by='gene',all.x=T)
SL_table<-merge(SL_table,SL,by='gene',all.x=T)
OS_table<-merge(OS_table,OS,by='gene',all.x=T)
ZM_table<-merge(ZM_table,ZM,by='gene',all.x=T)

AT_table$conserve[is.na(AT_table$conserve)]<-'1'
NB_table$conserve[is.na(NB_table$conserve)]<-'1'
SL_table$conserve[is.na(SL_table$conserve)]<-'1'
OS_table$conserve[is.na(OS_table$conserve)]<-'1'
ZM_table$conserve[is.na(ZM_table$conserve)]<-'1'

write.table(AT_table,'AT.tsv',col.names=T,row.names=T,sep='\t',quote=F)
write.table(NB_table,'NB.tsv',col.names=T,row.names=T,sep='\t',quote=F)
write.table(SL_table,'SL.tsv',col.names=T,row.names=T,sep='\t',quote=F)
write.table(OS_table,'OS.tsv',col.names=T,row.names=T,sep='\t',quote=F)
write.table(ZM_table,'ZM.tsv',col.names=T,row.names=T,sep='\t',quote=F)
