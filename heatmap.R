library(ggplot2)
library(reshape2)


myfunction <- function(input){
  file <- read.csv(input,sep=',',header = TRUE)
  filem <- melt(file)

  hp1=ggplot(filem, aes(variable, Protein, fill= value)) + 
    xlab("Samples")+
    ylab('Proteins')+
    geom_tile()+
    scale_fill_gradient(low="white", high="dark blue")+
    theme_bw()+
    theme(panel.border = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
    theme(axis.ticks = element_blank())
  hp1
}

myfunction('angiogenesis.csv')
myfunction('cytokine.csv')
