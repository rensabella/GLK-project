# Read the nucleootic diversity data group by group.
files<-dir('../../data/nucleotide_diversity_data')
plt_data<-NULL
for (file in files){
    group<-gsub('.txt','',unlist(strsplit(file,'_'))[2],fixed=T)
    pi<-read.table(paste0('../../data/nuleotide_diversity_data/',file))[,9]
    data<-data.frame(pi,group=rep(group,length(pi)))
    plt_data<-rbind(plt_data,data)
}
# Divide the data.
facet<-as.character(plt_data[,2])
facet[facet %in% c('ATAC','ChIP','genome')]<-'Binding'
facet[facet %in% c('X1','X2','X3','X4','X5')]<-'Conserved'
facet[facet %in% c('DEG','nonDEG','x54DEG','x54nonDEG')]<-'DEG'
plt_data$facet<-facet
plt_data$group<-factor(plt_data$group,levels=c('genome','ATAC','ChIP','X1','X2','X3','X4','X5','DEG','nonDEG','x54DEG','x54nonDEG'))

# Divide the data.
# ATAC vs ChIP
plt_genome_ATAC_CHIP<-plt_data[plt_data[,2] %in% c('ATAC','ChIP'),]
plt_genome_ATAC_CHIP$group<-factor(plt_genome_ATAC_CHIP$group,levels=c('ATAC','ChIP'))
# X1 vs X5
plt_conserve<-plt_data[plt_data[,2] %in% c('X1','X5'),]
plt_conserve$group<-factor(plt_conserve$group,levels=c('X1','X5'))
# DEG vs non-DEG
plt_DEG<-plt_data[plt_data[,2] %in% c('DEG','nonDEG'),]
plt_DEG$group<-factor(plt_DEG$group,levels=c('DEG','nonDEG'))
# DEG in conserve group vs non-DEG in conserve group
plt_x54_DEG<-plt_data[plt_data[,2] %in% c('x54DEG','x54nonDEG'),]
plt_x54_DEG$group<-factor(plt_x54_DEG $group,levels=c('x54DEG','x54nonDEG'))

# ATAC vs ChIP
x<-plt_genome_ATAC_CHIP[plt_genome_ATAC_CHIP[,2]=='ATAC',1]
y<-plt_genome_ATAC_CHIP[plt_genome_ATAC_CHIP[,2]=='ChIP',1]
pvalue<-paste0('pvalue==',ks.test(x, y)$p.value)
pdf('ATAC_CHIP.pdf',height=4,width=4)
ggplot(plt_genome_ATAC_CHIP, aes(x=pi, col=group)) + stat_ecdf(geom="smooth", se=F, size=1.2)+
theme_bw()+theme(
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
labs(x='Nucleotide diversity in ATAC and ChIP',y='Ecef (pi)')+scale_color_manual(values=c('#E69F00','#999999'))+
annotate("text",x=0.125,y=0.3,label=pvalue,parse=T)
dev.off()

# X1 vs X5
x<-plt_conserve[plt_conserve[,2]=='X1',1]
y<-plt_conserve[plt_conserve[,2]=='X5',1]
pvalue<-paste0('pvalue==',ks.test(x, y)$p.value)
pdf('X1_X5.pdf',height=4,width=4)
ggplot(plt_conserve, aes(x=pi, col=group)) + stat_ecdf(geom="smooth", se=F, size=1.2)+
theme_bw()+theme(
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
labs(x='Nucleotide diversity in X1 and X5',y='Ecef (pi)')+scale_color_manual(values=c('#E69F00','#999999'))+
annotate("text",x=0.03,y=0.3,label=pvalue,parse=T)
dev.off()

# DEG vs non-DEG
x<-plt_DEG[plt_DEG[,2]=='DEG',1]
y<-plt_DEG[plt_DEG[,2]=='nonDEG',1]
pvalue<-paste0('pvalue==',ks.test(x, y)$p.value)
pdf('DEG_ nonDEG.pdf',height=4,width=4)
ggplot(plt_DEG, aes(x=pi, col=group)) + stat_ecdf(geom="smooth", se=F, size=1.2)+
theme_bw()+theme(
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
labs(x='Nucleotide diversity in DEG and nonDEG',y='Ecef (pi)')+scale_color_manual(values=c('#E69F00','#999999'))+
annotate("text",x=0.065,y=0.3,label=pvalue,parse=T)
dev.off()

