# Read the duplication pattern data.
exp<-read.table('../../data/duplication_pattern_data/ZM_SG_EXP.tsv',header=T,sep='\t')
both_glk<-read.table('../../data/duplication_pattern_data/pair_54.txt')
sg1_glk<-read.table('../../data/duplication_pattern_data/pair_orpha_sg1_215')
sg2_glk<-read.table('../../data/duplication_pattern_data/pair_orpha_sg2_162')

both_glk<-as.matrix(both_glk)
sg1_glk<-as.matrix(sg1_glk)
sg2_glk<-as.matrix(sg2_glk)

# Calculate the correlation of paired genes.
both_glk_cor<-apply(both_glk,1,function(x){
cor(as.numeric(exp[exp[,1] %in% x[1],4:26]),as.numeric(exp[exp[,1] %in% x[2],4:26]))
})
both_glk_cor<-abs(na.omit(both_glk_cor))

sg1_glk_cor<-apply(sg1_glk,1,function(x){
cor(as.numeric(exp[exp[,1] %in% x[1],4:26]),as.numeric(exp[exp[,1] %in% x[2],4:26]))
})
sg1_glk_cor<-abs(na.omit(sg1_glk_cor))

sg2_glk_cor<-apply(sg2_glk,1,function(x){
cor(as.numeric(exp[exp[,1] %in% x[1],4:26]),as.numeric(exp[exp[,1] %in% x[2],4:26]))
})
sg2_glk_cor<-abs(na.omit(sg2_glk_cor))

# Rearrange the data to the shape that could be plotted.
both_plot_data<-data.frame(SG1=log2(exp[exp[,1] %in% both_glk[,1],15]),SG2=log2(exp[exp[,1] %in% both_glk[,2],15]))
both_plot_data<-both_plot_data[order(abs(both_plot_data[,2]-both_plot_data[,1]),decreasing=T),]
both_plot_data<-data.frame(group=1:nrow(both_plot_data),both_plot_data)

sg1_plot_data<-t(apply(sg1_glk,1,function(x){
c(log2(as.numeric(exp[exp[,1] %in% x[1],15])),log2(as.numeric(exp[exp[,1] %in% x[2],15])))}))
sg1_plot_data<-data.frame(sg1_plot_data)
colnames(sg1_plot_data)<-c('SG1','SG2')
sg1_plot_data<-sg1_plot_data[order(abs(sg1_plot_data[,2]-sg1_plot_data[,1]),decreasing=T),]
sg1_plot_data<-data.frame(group=1:nrow(sg1_plot_data),sg1_plot_data)

sg2_plot_data<-t(apply(sg2_glk,1,function(x){
c(log2(as.numeric(exp[exp[,1] %in% x[1],15])),log2(as.numeric(exp[exp[,1] %in% x[2],15])))}))
sg2_plot_data<-data.frame(sg2_plot_data)
colnames(sg2_plot_data)<-c('SG1','SG2')
sg2_plot_data<-sg2_plot_data[order(abs(sg2_plot_data[,2]-sg2_plot_data[,1]),decreasing=T),]
sg2_plot_data<-data.frame(group=1:nrow(sg2_plot_data),sg2_plot_data)

# Plot the correlation.
library(ggplot2)
library(see)
plot_data<-data.frame(correlation=c(both_glk_cor,sg1_glk_cor,sg2_glk_cor),group=c(rep('Both GLK',length(both_glk_cor)),rep('SG1 GLK',length(sg1_glk_cor)),rep('SG2 GLK',length(sg2_glk_cor))))

pdf('SG_raincloud.pdf')
ggplot(plot_data, aes(x = group,
                 y = correlation, 
                 fill = group)) +
  geom_violindot(alpha=0.5,color=NA,size=0.6) +
  theme_modern()+
  coord_flip()+
  scale_fill_social()
dev.off()
