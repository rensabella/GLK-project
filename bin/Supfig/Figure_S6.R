for (s1 in c('AT','NB','SL','OS','ZM')){
for (s2 in c('AT','NB','SL','OS','ZM')){
if (s1 != s2){
s1_data<-read.table(paste0('plot_2_species_cor/promoter_coverage_',s1,'_',s2,'_',s1))
s2_data<-read.table(paste0('plot_2_species_cor/promoter_coverage_',s1,'_',s2,'_',s2))
s1_s2_gene<-read.table(paste0('plot_2_species_cor/promoter_site_',s1,'_',s2,'_blast'))
s1_gene<-read.table(paste0(s1,'.tsv'),sep='\t',quote='')
s2_gene<-read.table(paste0(s2,'.tsv'),sep='\t',quote='')
group<-rep(paste0(s1,' target'),nrow(s1_s2_gene))
group[intersect(which(s1_s2_gene[,1] %in% s1_gene[,1]),which(s1_s2_gene[,2] %in% s2_gene[,1]))]<-'Both target'
group<-factor(group,levels=c('Both target',paste0(s1,' target')))
plot_data<-data.frame(s1=rowSums(s1_data[,4:5]),s2=rowSums(s2_data[,4:5]),group=group)

p<-ggplot(plot_data,aes(log2(s1),log2(s2),color=group))+geom_point(stroke=1.2,alpha=0.5,size=0.2)+geom_smooth(method='lm')+
theme_bw()+scale_color_manual(values=c('#E18727FF','#6F99ADFF'))+theme(panel.grid=element_blank())+
labs(x=paste0(s1,' log2 coverage'),y=paste0(s2,' log2 coverage'))
ggsave(paste0(s1,'_',s2,'_blast.pdf'),p,width=4,height=3)
}	
}
}