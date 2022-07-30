# Functional GO enrichment analysis and plot the top 5 results.
# The GO result files are got from agriGO v2.0.

library(ggplot2)
library(patchwork)

K0_AgriGO <- read.delim("../../data/GO_data/AT_result.txt", header = TRUE)
colnames(K0_AgriGO)
colnames(K0_AgriGO)[1] <- "term_ID"
colnames(K0_AgriGO)[3] <- "description"
K0_AgriGO_p005<-K0_AgriGO[K0_AgriGO$pvalue <=0.05,]
write.table(K0_AgriGO_p005,'AT_05.tsv',col.names=T,row.names=F,sep='\t',quote=F)

K0_REVIGO <- read.csv("../../data/GO_data/AT_Revigo.tsv", header = TRUE,sep='\t')
colnames(K0_REVIGO)[1]<-"term_ID" 
colnames(K0_REVIGO)[2]<-"description" 
K0_REVIGO_True <- K0_REVIGO[K0_REVIGO$Representative!='null',]

K0_GOenrich <- merge(x = K0_AgriGO_p005[,1:9], y = K0_REVIGO_True, by = c("term_ID","description"), all.y =TRUE)
K0_GOenrich$plotval <- -1 * log10(K0_GOenrich$pvalue)

K0_GOenrich$description <- factor(K0_GOenrich$description, levels = K0_GOenrich$description[order(K0_GOenrich$plotval)])
top5_one_K0 <- K0_GOenrich[order(K0_GOenrich$plotval,decreasing=T),][1:5,]

p<-ggplot(top5_one_K0, aes(x=description, y=plotval)) +
  geom_bar(stat="identity", width = 0.9) + coord_flip() +
  geom_text(aes(label=pvalue), vjust=0) +
  theme(axis.title = element_blank(),
        axis.text.y = element_text(face="bold", size=12),
        axis.text.x = element_blank(),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "white"),
        panel.background = element_blank()
        )+scale_y_discrete(breaks=c(0,25),expand = c(0,0))+geom_blank(aes(y=25))


