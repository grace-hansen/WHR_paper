#!/usr/bin/Rscript
library(optparse)
library(tidyverse)
library(data.table)
library(RCurl)
arguments <- parse_args(OptionParser(usage = "make_sex_scatterplot.R <trait> <tissue> <type> <data_source>",option_list=list()),
                        positional_arguments = 4)
opt<-arguments$opt
trait<-arguments$args[1]
tissue<-arguments$args[2]
type<-arguments$args[3]
data_source<-arguments$args[4]
setwd("~/medusa/papers/TWAS/TWAS/sex_scatter/")
pony_colors<-fread("~/medusa/papers/TWAS/pony_palette")

######### Author: Grace Hansen #########
#This script plots effect sizes from male and female TWAS in a scatterplot/..

########## For color manipulation of graph ###############
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

darken <- function(color, factor=1.2){
  col <- col2rgb(color)
  col <- col/factor
  col <- rgb(t(col), maxColorValue=255)
  col
}

########## Gather data: get list of all genes, add column denoting significance
if (type=="expression") {
  F_res<-fread(paste("~/midway/",type,"/",trait,"_F/",tissue,"/results/",data_source,".all.dat",sep=''))
  M_res<-fread(paste("~/midway/",type,"/",trait,"_M/",tissue,"/results/",data_source,".all.dat",sep=''))
} else if (type=="splicing") {
  F_res<-fread(paste("~/midway/",type,"/",trait,"_F/",tissue,"/results/",data_source,".all.dat.gene",sep='')) %>%
    filter(.,GENE!="<NA>") %>%
    group_by(GENE) %>% top_n(-1, TWAS.P)
  M_res<-fread(paste("~/midway/",type,"/",trait,"_M",tissue,"/results/",data_source,".all.dat.gene",sep='')) %>%
    filter(.,GENE!="<NA>") %>%
    group_by(GENE) %>% top_n(-1, TWAS.P)
  F_res$ID=F_res$GENE
  M_res$ID=M_res$GENE
} else if (type=="both") {
  
}

#Incorporate significance
sig_genes_F=scan(paste("~/midway/",type,"/",trait,"_F/",tissue,"/results/posthoc/sig_genes",sep=''),what="character",sep='\n')
sig_genes_M=scan(paste("~/midway/",type,"/",trait,"_M/",tissue,"/results/posthoc/sig_genes",sep=''),what="character",sep='\n')
for (i in 1:nrow(F_res)) {
  if (F_res$ID[i] %in% sig_genes_F) {
    F_res$sig[i]<-"S"
  } else {
    F_res$sig[i]<-"NS"
  }
}
for (i in 1:nrow(M_res)) {
  if (M_res$ID[i] %in% sig_genes_M) {
    M_res$sig[i]<-"S"
  } else {
    M_res$sig[i]<-"NS"
  }
}

############### Merge male and female data ##################
dat<-merge(F_res,M_res,by="ID",suffixes=c("_F","_M"))
dat$sig=NA
for (i in 1:nrow(dat)) {
  if (dat$sig_F[i]=="S" && dat$sig_M[i]=="S") {
    dat$sig[i]="both"
  } else if (dat$sig_F[i]=="S" && dat$sig_M[i]=="NS") {
    dat$sig[i]="female"
  } else if (dat$sig_F[i]=="NS" && dat$sig_M[i]=="S") {
      dat$sig[i]="male"
  } else if (dat$sig_F[i]=="NS" && dat$sig_M[i]=="NS") {
    dat$sig[i]="neither" }
}

############ Plot data ###############
F_abs_limit=min(abs(head(dat[order(dat$TWAS.Z_F),],5)$TWAS.Z_F)) #Get top 10 female results for labelling
M_abs_limit=min(abs(head(dat[order(dat$TWAS.Z_M),],5)$TWAS.Z_M)) #Get top 10 male results for labelling
G<-ggplot(dat,aes(x=abs(dat$TWAS.Z_F),y=abs(dat$TWAS.Z_M),color=dat$sig))+
  geom_point(size=2.5)+
  geom_text(data=dat[(abs(dat$TWAS.Z_F) >= F_abs_limit | abs(dat$TWAS.Z_M) >= M_abs_limit),],
            aes(abs(TWAS.Z_F),abs(TWAS.Z_M),label=ID,color=sig),
            position=position_nudge(y=0.4),
            size=3.5,fontface="bold.italic")+
  labs(color="Significance")+
  scale_color_manual(values = c(rgb(pony_colors[11,1:3]),darken(rgb(pony_colors[10,1:3])),darken(rgb(pony_colors[14,1:3])),as.character(add.alpha("gray60",alpha=0.2))))+
  scale_y_continuous(name="Male absolute value TWAS Z scores", limits=c(0,12))+
  scale_x_continuous(name="Female absolute value TWAS Z scores",limits=c(0,24))+
  theme_minimal()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14))
G


pdf(paste("~/medusa/papers/TWAS/TWAS/sex_scatter/",trait,"_",tissue,"_",type,"_sex_scatterplot.pdf",sep=''),width=8,height=3.25)
G
dev.off()

