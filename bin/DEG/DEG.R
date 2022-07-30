library(DESeq2)
# Find out DEG with package DESeq2.
# The gene count results are calculate with hisat-count.
count_files <- dir()[grep('SRR.*.exon.*count.txt',dir())]
deg_table <- rep(0,35176)
for (file in count_files){

    data <- read.table(file)
    deg_table <- data.frame(deg_table,data[,2])

}
deg_table <- deg_table[,-1]
colnames(deg_table) <- c('mut_1','mut_2','mut_3','wild_1','wild_2','wild_3')
rownames(deg_table) <- data[,1]
condition  <-  factor(c(rep("mut",3), rep("wild",3)))    
coldata  <-  data.frame(row.names = colnames(deg_table), condition)

# DESeq2
dds  <-  DESeqDataSetFromMatrix(countData=deg_table, colData=coldata, design=~condition)
dds_uni <- DESeq(dds)    
resultsNames(dds_uni)    
dds_uni$condition        
res <- results(dds_uni)  
summary(res)        
table(res$padj < 0.05)        
res <- res[order(res$padj),]  
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds_uni, normalized=TRUE)),by="row.names",sort=FALSE)
deg_genes_result <- resdata[which(resdata$padj<0.05),]

# Combine the FPKM value calculated by cufflinks
fpkm_files <- dir()[grep('SRR.*.cufflinks$',dir())]
for (f in fpkm_files){
    fpkm <- read.table(paste0(f,'/isoforms.fpkm_tracking'),header=T)
    resdata <- merge(resdata,fpkm[,c(1,10)],by.x='Row.names',by.y='tracking_id')
}
colnames(resdata)[14:19] <- paste0('FPKM_',colnames(resdata)[8:13])
multi_iso_genes <- unique(substr(resdata[which(duplicated(substr(resdata$Row.names,1,9))),1],1,9))
for (gene in multi_iso_genes){
    gene_site <- grep(gene,resdata$Row.names)
    delete_isoform_site <- gene_site[order(resdata[gene_site,]$FPKM_mut_1,decreasing=T)][-1]
    resdata <- resdata[-delete_isoform_site,]
}
resdata$Row.names <- substr(resdata$Row.names,1,9)
write.table(resdata,file = "new_deg_result_fpkm.tsv",col.names=T,row.names=F,sep='\t',quote=F)
