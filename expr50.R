expr50 <- read.csv('sample_exprSet_14450_nrALL.txt',sep = '\t', header = TRUE)

expr50Set <- expr50[,c("ID_REF","GSM1149950","GSM1149952","GSM1149985","GSM1149993","GSM1149996","GSM1150010","GSM1150056","GSM1150067","GSM1150089","GSM1150100","GSM1150119","GSM1150122","GSM1150169","GSM1150173","GSM1150177","GSM1150179","GSM1150183","GSM1150191","GSM1150196","GSM1150197","GSM1150198","GSM1150207","GSM1150208","GSM1150209","GSM1150221","GSM1150229","GSM1150231","GSM1150232","GSM1150233","GSM1150244","GSM1150245","GSM1150247","GSM1150256","GSM1150260","GSM1150264","GSM1150267","GSM1150270","GSM1150277","GSM1150280","GSM1150281","GSM1150286","GSM1150302","GSM1150306","GSM1150313","GSM1150354","GSM1150355","GSM1150363","GSM1150365","GSM1150367","GSM1150370","GSM1150372","GSM1150375","GSM1150377","GSM1150378","GSM1150379","GSM1150380","GSM1150382","GSM1150384","GSM1150385","GSM1150387","GSM1150390","GSM1150391","GSM1150392","GSM1150395","GSM1150398","GSM1150399","GSM1150400","GSM1150401","GSM1150404","GSM1150410","GSM1150412","GSM1150415","GSM1150416","GSM1150418","GSM1150421","GSM1150424","GSM1150425","GSM1150426","GSM1150428","GSM1150429","GSM1150431","GSM1150432","GSM1150436","GSM1150441","GSM1150442","GSM1150445","GSM1150446","GSM1150448","GSM1150449","GSM1150451","GSM1150454","GSM1150460","GSM1150461","GSM1150466","GSM1150467","GSM1150469","GSM1150470","GSM1150475","GSM1150480","GSM1150484","GSM1150487","GSM1150488","GSM1150489","GSM1150490","GSM1150495","GSM1150497","GSM1150503","GSM1150504","GSM1150510","GSM1150513","GSM1150519","GSM1150525","GSM1150530","GSM1150532","GSM1150533","GSM1150534","GSM1150535","GSM1150536","GSM1150541","GSM1150542","GSM1150545","GSM1150546","GSM1149948","GSM1149951","GSM1149958","GSM1149962","GSM1149963","GSM1149966","GSM1149971","GSM1149975","GSM1149981","GSM1149987","GSM1149988","GSM1149992","GSM1149998","GSM1150001","GSM1150012","GSM1150014","GSM1150017","GSM1150030","GSM1150031","GSM1150032","GSM1150037","GSM1150040","GSM1150041","GSM1150053","GSM1150055","GSM1150070","GSM1150075","GSM1150078","GSM1150080","GSM1150081","GSM1150085","GSM1150090","GSM1150143","GSM1150144","GSM1150146","GSM1150148","GSM1150149","GSM1150150","GSM1150151","GSM1150152","GSM1150153","GSM1150154","GSM1150156","GSM1150157","GSM1150158","GSM1150159","GSM1150160","GSM1150161","GSM1150165","GSM1150167","GSM1150168","GSM1150171","GSM1150174","GSM1150175","GSM1150176","GSM1150373","GSM1150376","GSM1150381","GSM1150394","GSM1150402","GSM1150403","GSM1150405","GSM1150413","GSM1150417","GSM1150423","GSM1150433","GSM1150435","GSM1150438","GSM1150443","GSM1150444","GSM1150447","GSM1150452","GSM1150457","GSM1150459","GSM1150463","GSM1150468","GSM1150473","GSM1150478","GSM1150491","GSM1150498","GSM1150501","GSM1150506","GSM1150509","GSM1150515","GSM1150516","GSM1150517","GSM1150526","GSM1150527","GSM1150531","GSM1150537","GSM1150543")]

expr50Setsub <- as.matrix(expr50Set[,2:214])
rownames(expr50Setsub) <- expr50Set$ID_REF

par(cex = 0.7)
n.expr50Setsub=ncol(expr50Setsub)
if(n.expr50Setsub>214) par(cex = 0.5)
cols <- rainbow(n.expr50Setsub*1.2)
boxplot(expr50Setsub, col = cols,main="expression value",las=2)

group50 <- read.csv('group_list.txt',sep = '\t',header = FALSE)
design50 <- model.matrix(~0+factor(group50$V2))
colnames(design50) <- levels(factor(group50$V2))
rownames(design50) <- group50$V1

contrast50 <- makeContrasts(IPF-CTRL,levels = design50)

fit <- lmFit(expr50Setsub,design50)
fit2 <- contrasts.fit(fit, contrast50)
fit2 <- eBayes(fit2) 
tempOutput <- topTable(fit2, coef=1, n=Inf)
nrDEG50 <- na.omit(tempOutput)

library(ggplot2)
library(ggrepel)
nrDEG50$threshold <- factor(ifelse(nrDEG50$adj.P.Val < 0.05 & abs(nrDEG50$logFC) >= 1, 
                                 ifelse(nrDEG50$logFC>= 1 ,'Up','Down'),'NoSignificant'),
                          levels=c('Up','Down','NoSignificant'))

nrDEG50$target <- ''
tar <- c(1924,7099,808,9,91,98,89,7,457)
#View(tar)
nrDEG50$target[tar] <- rownames(nrDEG50)[tar]
#View(nrDEG$target)
#View(nrDEG50)
#---------------------------Volcano2----------------------------------------
p2 <- ggplot(nrDEG50,aes(x=logFC,y=-log10(adj.P.Val), label = target,color=threshold))+
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
    data          = nrDEG50,
    nudge_y       = 30 - nrDEG50$adj.P.Val,
    segment.size  = 0.5,
    segment.color = "purple",
    direction     = "x",
    show.legend = FALSE
  ) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p2
