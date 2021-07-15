#!/usr/bin/Rscript
args=commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("USAGE: plot_top_genes.R <gene or intron list> <trait> <tissue> <type> <data_source> <number of genes> \n", call.=FALSE)
}
library(data.table)
library(tidyverse)

trait=args[1]
tissue=args[2]
type=args[3]
data_source=args[4]
topN=as.numeric(args[5])
pony_colors<-fread("~/medusa/papers/TWAS/pony_palette")
setwd("~/medusa/papers/TWAS/TWAS/top_genes/")

######### Author: Grace Hansen #########
#This script plots the top 10 genes from a TWAS.

###### Adjust graph colors ###########
darken <- function(color, factor=1.2){
  col <- col2rgb(color)
  col <- col/factor
  col <- rgb(t(col), maxColorValue=255)
  col
}
######################################

############ Gather data ###############
if (type=="expression") {
  res<-fread(paste("~/midway/",type,"/",trait,"/",tissue,"/results/",data_source,".all.dat.top",sep=''))
} else if (type=="splicing") {
  res<-fread(paste("~/midway/",type,"/",trait,"/",tissue,"/results/",data_source,".all.dat.top.gene",sep='')) 
  res$ID=res$GENE
}
pthresh=0.05/nrow(res)
res<-res %>%
  filter(.,ID!="<NA>") %>%
  group_by(ID) %>% filter(., rank(TWAS.P, ties.method="first")==1)

res<-res[1:topN,]
res$TWAS.Z<-as.numeric(res$TWAS.Z)
res$ID<-factor(res$ID,levels=as.character(rev(res$ID)))

########## Plot #############
if (trait=="obesity") {
  pony=darken(rgb(pony_colors[3,1:3]))
} else if (trait=="WHR_F") {
  pony=darken(rgb(pony_colors[10,1:3]))
} else if (trait=="WHR_M") {
  pony=rgb(pony_colors[14,1:3])
}
pdf(paste("/home/grace/medusa/papers/TWAS/TWAS/top_genes/",trait,"_",tissue,"_",type,"_top_genes.pdf",sep=''),width=4,height=topN/5+1)

ggplot(data=res,aes(x=ID,y=abs(TWAS.Z)))+
  geom_point(stat="identity",size=5,color=pony)+
  geom_hline(yintercept=0,lwd=1)+
  scale_y_continuous(name="Absolute TWAS Z score",limits=c(0,30))+
  theme_minimal()+
  theme(axis.text=element_text(size=19),
        axis.text.y=element_text(face="italic"),
        axis.title.y=element_blank(),
        axis.title.x=element_text(size=19))+
  coord_flip()

dev.off()
