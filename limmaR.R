#--------------------------Input-------------------------------------------
exp2 <- read.table('limma_sample-matrix-IPF-CTRL-6480-sort.txt',sep = '\t',header = TRUE)
data <- exp2[,2:56]
data
#---------------Expression matrix---------------------
data <- as.matrix(data)
rownames(data) <- exp2$X
#-----------------------------------------------------

#-------------------------QC-------------------------
#View(data)
par(cex = 0.7)
n.data=ncol(data)
if(n.data>200) par(cex = 0.5)
cols <- rainbow(n.data*1.2)
boxplot(data, col = cols,main="expression value",las=2)
#---------------------------------------------------------

#------------------------Design matrix-------------------------------
group_list <- read.csv('group_list.txt',sep = '\t',header = FALSE)
design <- model.matrix(~0+factor(group_list$V2))
colnames(design) <- levels(factor(group_list$V2))
rownames(design) <- group_list$V1
#--------------------------------------------------------------------

#---------------------Contrast Matrix--------------------------------
contrast <- makeContrasts(IPF-Control,levels = design)
#--------------------------------------------------------------------

#--------------------------limma fitting-----------------------------
fit <- lmFit(data,design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2) 
tempOutput <- topTable(fit2, coef=1, n=Inf)
nrDEG <- na.omit(tempOutput) 
#head(nrDEG)
#----------------------------output---------------------------------
write.csv(nrDEG,"GPL6480_IPFvsContrl_diffExpr.csv",quote = F)

#-------------------------Volcanol plot------------------------------
library(ggplot2)
library(ggrepel)
#-----------------------------Vocano1----------------------------
nrDEG$color <- ifelse(nrDEG$adj.P.Val<0.05 & abs(nrDEG$logFC)>= 2,
                    ifelse(nrDEG$logFC > 2,'red','blue'),'gray')
color <- c(red = "red",gray = "gray",blue = "blue")

p1 <- ggplot(nrDEG, aes(x=logFC, y=-log10(adj.P.Val),col = color))+
  geom_point()+
  theme_bw() +
  scale_color_manual(values = color) +
  labs(x="log2 (fold change)",y="-log10 (q-value)") +
  geom_hline(yintercept = -log10(0.05), lty=4,col="grey",lwd=0.6) +
  geom_vline(xintercept = c(-2, 2), lty=4,col="grey",lwd=0.6) +
  theme(legend.position = "none",
        panel.grid=element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))
p1
#----------------------Volcano2_label-------------------------------------------------------
nrDEG$threshold <- factor(ifelse(nrDEG$adj.P.Val < 0.05 & abs(nrDEG$logFC) >= 1, 
                                 ifelse(nrDEG$logFC>= 1 ,'Up','Down'),'NoSignificant'),
                          levels=c('Up','Down','NoSignificant'))
nrDEG$target <- ''
tar <- c(4,284,1048,1443,8214)
#View(tar)
nrDEG$target[tar] <- rownames(nrDEG)[tar]
#View(nrDEG$target)
#View(nrDEG)
#---------------------------Volcano2----------------------------------------
p2 <- ggplot(nrDEG,aes(x=logFC,y=-log10(adj.P.Val), label = target,color=threshold))+
  geom_point()+
  scale_color_manual(values=c("#DC143C","#00008B","#808080"))+#确定点的颜色
  theme_bw()+#修改图片背景
  theme(
    legend.title = element_blank()#不显示图例标题
  )+
  ylab('-log10 (p-adj)')+#修改y轴名称
  xlab('log2 (FoldChange)')+#修改x轴名称
  geom_vline(xintercept=c(-1,1),lty=3,col="black",lwd=0.5) +#添加横线|FoldChange|>2
  geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=0.5)+#添加竖线padj<0.05
  geom_text_repel(
    data          = nrDEG,
    nudge_y       = 30 - nrDEG$adj.P.Val,
    segment.size  = 0.5,
    segment.color = "purple",
    direction     = "x",
    show.legend = FALSE
  ) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p2

