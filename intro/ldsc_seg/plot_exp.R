#!/usr/bin/Rscript
library(optparse)
library(tidyverse)
library(data.table)
arguments <- parse_args(OptionParser(usage = "plot_exp.R <trait>",option_list=list()),
                        positional_arguments = 1)
opt<-arguments$opt
trait<-arguments$args[1]

setwd(paste("/home/grace/midway/ldsc_seg",sep=""))
pony_colors<-fread("~/medusa/papers/TWAS/pony_palette")

########### Author: Grace Hansen ###############
# This script plots enrichment for specifically-espressed genes from ldsc-seg.

# Load data
dat<-fread(paste("output/",trait,"/GTEx_exp/all_tissues.results",sep=''))
groups<-fread(paste("output/",trait,"/GTEx_exp/tissues_groups",sep=''))
dat<-merge(groups,dat,by="Tissue")
dat<-dat %>% arrange(Group,-log10(Enrichment_p))
dat$Tissue<-factor(dat$Tissue,levels=unique(dat$Tissue))

set.seed(29)
pony_palette=slice(pony_colors,sample(1:length(unique(dat$Group))))

#Plot mean of each group
plot_dat<-dat %>% group_by(Group) %>% summarize(p=mean(Enrichment_p))
plot_dat$Group<-as.factor(plot_dat$Group)
G<-ggplot(plot_dat,aes(x=Group,y=-log10(p)))+
  geom_bar(aes(fill=Group),position="dodge",stat="identity")+
  theme_minimal()+
  scale_y_continuous("mean -log10 Enrichment p-value")+
  scale_fill_manual(values=rgb(pony_palette[,1:3]))+
  theme(axis.title=element_blank(),
        axis.text=element_text(size=22),
        legend.position = "none")+
  scale_x_discrete(breaks = NULL,limits=rev(levels(plot_dat$Group)))

pdf(paste("output/",trait,"/GTEx_exp/GTEx_ldscseg_exp_groupmean.pdf",sep=''),width=5,height=3)
G
dev.off()


#Plot each member of each group
G<-ggplot(dat,aes(x=Tissue,y=-log10(Enrichment_p)))+
  geom_bar(aes(fill=Group),position="dodge",stat="identity")+
  theme_minimal()+
  scale_y_continuous("-log10 Enrichment p-value")+
  scale_fill_manual(values=rgb(pony_palette[,1:3]))+
  theme(axis.title.y=element_text(size=15),
        axis.title.x=element_blank(),
        axis.text.x=element_text(size=10,angle=-90),
        legend.position="none")


pdf(paste("output/",trait,"/GTEx_exp/GTEx_ldscseg_exp.pdf",sep=''),width=7.5,height=4.5)
G
dev.off()
