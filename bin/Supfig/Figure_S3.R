## Arabidopsis
library(Rmisc)
library(ggpubr)
library(ggplot2)

data <- read.table("Readscount_AtGLK_500bp.txt", header = TRUE)
countdata <- data[,c(1,2,3,4,8,12,16)]
colnames(countdata) <- c("Chr","start","end","AtGLK1_rep1","AtGLK1_rep2","AtGLK2_rep1","AtGLK2_rep2")

GLK1_1_sum <- sum(countdata$AtGLK1_rep1)
GLK1_2_sum <- sum(countdata$AtGLK1_rep2)
GLK2_1_sum <- sum(countdata$AtGLK2_rep1)
GLK2_2_sum <- sum(countdata$AtGLK2_rep2)

countdata$GLK1_rep1_CPM <- countdata$ZmGLK1_rep1/GLK1_1_sum*10000
countdata$GLK1_rep2_CPM <- countdata$ZmGLK1_rep2/GLK1_2_sum*10000
countdata$GLK2_rep1_CPM <- countdata$ZmGLK2_rep1/GLK2_1_sum*10000
countdata$GLK2_rep2_CPM <- countdata$ZmGLK2_rep2/GLK2_2_sum*10000

countdata$GLK1_CPM_mean <- (countdata$GLK1_rep1_CPM + countdata$GLK1_rep2_CPM)/2
countdata$GLK2_CPM_mean <- (countdata$GLK2_rep1_CPM + countdata$GLK2_rep2_CPM)/2
countdata$GLK1_CPM_log2 <- log2(countdata$GLK1_CPM_mean+1)
countdata$GLK2_CPM_log2 <- log2(countdata$GLK2_CPM_mean+1)

At_CPM_log <- ggplot(data = countdata, aes(x=GLK1_CPM_log2,y=GLK2_CPM_log2))+geom_point(size=0.5) 
    + geom_smooth(method=lm) + theme_bw() + theme(panel.grid=element_blank()) 
    + labs(x ='AtGLK1 log2CPM',y="AtGLK2 log2CPM",title = "AtGLK1/2 union peaks")
