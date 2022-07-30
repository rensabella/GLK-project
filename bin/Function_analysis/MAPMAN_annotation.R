# MAPMAN
# The MAPMAN analysis are performed with protein sequences.
AT<-read.table('../../data/mapman_data/Arabidopsis.results.txt',sep='\t',header=T)
NB<-read.table('../../data/mapman_data/Tobacco.results.txt',sep='\t',header=T)
SL<-read.table('../../data/mapman_data/Tomato.results.txt',sep='\t',header=T)
OS<-read.table('../../data/mapman_data/Rice.results.txt',sep='\t',header=T)
ZM<-read.table('../../data/mapman_data/Maize.results.txt',sep='\t',header=T)

# Uniform gene IDs.
AT[,3]<-gsub('g','G',gsub('at','AT',AT[,3]))
NB[,3]<-gsub('niben','Niben',gsub('chr','Chr',NB[,3]))
SL[,3]<-gsub('soly','Soly',SL[,3])
OS[,3]<-gsub('loc_o','LOC_O',OS[,3])
ZM[,3]<-gsub('zm','Zm',ZM[,3])
AT[,3]<-as.character(sapply(AT[,3],function(x) unlist(strsplit(x,'.',fixed=T))[1]))
NB[,3]<-as.character(sapply(NB[,3],function(x) unlist(strsplit(x,'.',fixed=T))[1]))
SL[,3]<-as.character(sapply(SL[,3],function(x) unlist(strsplit(x,'.',fixed=T))[1]))
OS[,3]<-as.character(sapply(OS[,3],function(x) unlist(strsplit(x,'.',fixed=T))[1]))
ZM[,3]<-as.character(sapply(ZM[,3],function(x) unlist(strsplit(x,'_',fixed=T))[1]))

# The GLK target IDs inforamtion is in Supplementary Table S3.
AT_glk_target<-read.table('AT.tsv',sep='\t',quote='')[,1]
NB_glk_target<-read.table('NB.tsv',sep='\t',quote='')[,1]
SL_glk_target<-read.table('SL.tsv',sep='\t',quote='')[,1]
OS_glk_target<-read.table('OS.tsv',sep='\t',quote='')[,1]
ZM_glk_target<-read.table('ZM.tsv',sep='\t',quote='')[,1]

# Analysis with hypergeometric enrichment.
AT_result<-NULL
for (i in c(1:28,30)){
	term=AT[AT[,1]==i,2]
	term_genes=AT[grep(term,AT[,2]),3]
	k<-length(which(unique(term_genes) %in% AT_glk_target)) -1 
	N<-nrow(AT)
	M<-length(term_genes)
	n<-length(AT_glk_target)
    pvalue<-phyper(k, M, N-M, n, lower.tail=FALSE)
	AT_result<-rbind(AT_result,c(term, pvalue,k))
}

NB_result<-NULL
for (i in c(1:28,30)){
	term=NB[NB[,1]==i,2]
	term_genes=NB[grep(term,AT[,2]),3]
	k<-length(which(unique(term_genes) %in% NB_glk_target)) -1 
	N<-nrow(NB)
	M<-length(term_genes)
	n<-length(NB_glk_target)
    pvalue<-phyper(k, M, N-M, n, lower.tail=FALSE)
	NB_result<-rbind(NB_result,c(term, pvalue,k))
}

SL_result<-NULL
for (i in c(1:28,30)){
	term=SL[SL[,1]==i,2]
	term_genes=SL[grep(term,SL[,2]),3]
	k<-length(which(unique(term_genes) %in% SL_glk_target)) -1 
	N<-nrow(SL)
	M<-length(term_genes)
	n<-length(SL_glk_target)
    pvalue<-phyper(k, M, N-M, n, lower.tail=FALSE)
	SL_result<-rbind(SL_result,c(term, pvalue,k))
}

OS_result<-NULL
for (i in c(1:28,30)){
	term=OS[OS[,1]==i,2]
	term_genes=OS[grep(term,OS[,2]),3]
	k<-length(which(unique(term_genes) %in% OS_glk_target)) -1 
	N<-nrow(OS)
	M<-length(term_genes)
	n<-length(OS_glk_target)
    pvalue<-phyper(k, M, N-M, n, lower.tail=FALSE)
	OS_result<-rbind(OS_result,c(term, pvalue,k))
}

ZM_result<-NULL
for (i in c(1:28,30)){
	term=ZM[ZM[,1]==i,2]
	term_genes=ZM[grep(term,ZM[,2]),3]
	k<-length(which(unique(term_genes) %in% ZM_glk_target)) -1 
	N<-nrow(ZM)
	M<-length(term_genes)
	n<-length(ZM_glk_target)
    pvalue<-phyper(k, M, N-M, n, lower.tail=FALSE)
	ZM_result<-rbind(ZM_result,c(term, pvalue,k))
}

# Select the terms that term gene number >10.
select_n<-which(((AT_result[,3]>10)+(NB_result[,3]>10)+(SL_result[,3]>10)+(OS_result[,3]>10)+(ZM_result[,3]>10))==5)
AT_result<-AT_result[select_n,]
NB_result<-NB_result[select_n,]
SL_result<-SL_result[select_n,]
OS_result<-OS_result[select_n,]
ZM_result<-ZM_result[select_n,]

# Plot the p-value of the enrichment result.
library('ggplot2')
library('RColorBrewer')

mat<-rbind(AT_result,NB_result,SL_result,OS_result,ZM_result)
mat<-data.frame(GO= mat[,1],p_value=mat[,2],Number=mat[,3], 
Species=rep(c('Atha','Nben','Slyc','Osat','Zmay'),each=14))
mat$Species <- factor(mat$Species,levels=c('Atha','Nben','Slyc','Osat','Zmay'))
mat$GO <- factor(mat$GO,levels=unique(as.vector(mat[order(mat$p_value),]$GO)))
class(mat$p_value)<-'numeric'
class(mat$Number)<-'numeric'

myPalette <- colorRampPalette(rev(brewer.pal(11, "RdGy")))
pdf("Figure3c.pdf",height=6,width=5)
ggplot(mat,aes(Species,GO,size=Number,fill=-log10(p_value)))+
  geom_point(shape=21)+
  theme_bw()+
  scale_size(range=c(1,10))+
  scale_fill_gradientn(colours=myPalette(11))
dev.off()



library(pheatmap)
library(viridis)

# The GLK target IDs and conservation inforamtion is in Supplementary Table S3.
AT_glk_target<-read.table('AT.tsv',sep='\t',quote='')[,c('gene','conserve')]
NB_glk_target<-read.table('NB.tsv',sep='\t',quote='')[,c('gene','conserve')]
SL_glk_target<-read.table('SL.tsv',sep='\t',quote='')[,c('gene','conserve')]
OS_glk_target<-read.table('OS.tsv',sep='\t',quote='')[,c('gene','conserve')]
ZM_glk_target<-read.table('ZM.tsv',sep='\t',quote='')[,c('gene','conserve')]
terms<-unique(as.vector(mat[order(mat$p_value),]$GO))[14:1]

# Conservation Score calculation.
AT_conserve_score<-NULL
for (term in terms){
	term_genes=AT[grep(term,AT[,2]),3]
	GLK_term_genes = AT_glk_target[AT_glk_target[,1] %in% term_genes,]
    GLK_term_gene_score  = sum(GLK_term_genes[,2]/5)/length(GLK_term_genes[,2])
	AT_conserve_score <-rbind(AT_conserve_score,c(term, GLK_term_gene_score))
}
NB_conserve_score<-NULL
for (term in terms){
	term_genes=NB[grep(term,NB[,2]),3]
	GLK_term_genes = NB_glk_target[NB_glk_target[,1] %in% term_genes,]
    GLK_term_gene_score  = sum(GLK_term_genes[,2]/5)/length(GLK_term_genes[,2])
	NB_conserve_score <-rbind(NB_conserve_score,c(term, GLK_term_gene_score))
}
SL_conserve_score<-NULL
for (term in terms){
	term_genes=SL[grep(term,SL[,2]),3]
	GLK_term_genes = SL_glk_target[SL_glk_target[,1] %in% term_genes,]
    GLK_term_gene_score  = sum(GLK_term_genes[,2]/5)/length(GLK_term_genes[,2])
	SL_conserve_score <-rbind(SL_conserve_score,c(term, GLK_term_gene_score))
}
OS_conserve_score<-NULL
for (term in terms){
	term_genes=OS[grep(term,OS[,2]),3]
	GLK_term_genes = OS_glk_target[OS_glk_target[,1] %in% term_genes,]
    GLK_term_gene_score  = sum(GLK_term_genes[,2]/5)/length(GLK_term_genes[,2])
	OS_conserve_score <-rbind(OS_conserve_score,c(term, GLK_term_gene_score))
}
ZM_conserve_score<-NULL
for (term in terms){
	term_genes=ZM[grep(term,ZM[,2]),3]
	GLK_term_genes = ZM_glk_target[ZM_glk_target[,1] %in% term_genes,]
    GLK_term_gene_score  = sum(GLK_term_genes[,2]/5)/length(GLK_term_genes[,2])
	ZM_conserve_score <-rbind(ZM_conserve_score,c(term, GLK_term_gene_score))
}

# Plot the conservation score in GLK targets.
plt_data<-cbind(AT_conserve_score[,2], NB_conserve_score[,2], SL_conserve_score[,2], OS_conserve_score[,2], ZM_conserve_score[,2])
colnames(plt_data)<-c('Atha','Nben','Slyc','Osat','Zmay')
mode(plt_data)<-'numeric'

pdf("Fig3d.pdf",height=9,width=3)
pheatmap(plt_data,cluster_rows=F,cluster_cols=F,cellwidth=20,cellheight=20)+ scale_fill_distiller(palette = "Spectral")
dev.off()


