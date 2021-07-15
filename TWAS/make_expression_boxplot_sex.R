#!/usr/bin/Rscript
library(optparse)
library(tidyverse)
library(data.table)
library(gridExtra)

trait<-"WHR"
tissue<-"adipose_subcutaneous"
data_source<-"GTEx_v8"
gene<-"SNX10"
rsid<-"rs1534696"
setwd(paste("/home/grace/midway/expression/",trait,"/",tissue,"/results/posthoc",sep=""))
chrom=system(paste("zcat ../../",tissue,"_",data_source,"_exp.txt.gz | awk '{ if ($3 == \"",gene,"\") print $1}'",sep=''),intern=TRUE)
pony_colors<-fread("~/medusa/papers/TWAS/pony_palette")

#Make expression phenotpe file
if (data_source=="GTEx_v8") {
  tpms<-fread(cmd="zcat ~/midway/expression/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz | tail -n +3")
  samples<-scan(paste("~/midway/expression/",trait,"/",tissue,"/",tissue,"_sample_IDs",sep=''),what="character",sep="\n")
  colnames<-c("Description",samples)
  tpms<-tpms[,..colnames]
  tpms<-tpms[tpms$Description==gene]
  tpms$Description<-NULL
  subj_IDs<-paste("GTEX-",sapply(strsplit(names(tpms),'-'),'[[',2),sep='')
  names(tpms)<-subj_IDs
  write.table(tpms,"gene_pheno",row.names = FALSE,quote=FALSE)
  gc()
} else {
  cmd=paste("zcat ../../",tissue,"_",data_source,"_exp_1.txt.gz | head -1 > gene_pheno",sep='')
  system(cmd)
  cmd=paste("zcat ../../",tissue,"_",data_source,"_exp_1.txt.gz | awk '{ if ($3 == \"",gene,"\") print $0}' >> gene_pheno",sep='')
  system(cmd)
}
gene_dat<-fread("gene_pheno")
gene_dat<-as.data.frame(t(gene_dat),stringsAsFactors = FALSE)
rownames<-rownames(gene_dat)[4:length(rownames(gene_dat))]
gene_dat<-as.data.frame(gene_dat[4:nrow(gene_dat),])
colnames(gene_dat)<-"Expression"
IDs<-paste("GTEX-",sapply(strsplit(rownames,'-'),'[[',2),sep='')
gene_dat$P1<-IDs
gene_dat$P2<-IDs
gene_dat<-gene_dat %>% select(P1,P2,everything())
write.table(gene_dat,"gene.pheno",row.names = FALSE,col.names = FALSE,quote=FALSE)

#Add rsid genotypes to expression file from plink
cmd=paste("awk '{ if ($2 == \"",rsid,"\") print $4}' ~/midway/genos/",data_source,"/plink/chr",chrom,"_rsids.bim",sep="")
rsid_loc=as.numeric(system(cmd,intern=TRUE))
start=rsid_loc-5
stop=rsid_loc+5
cmd<-paste("plink --bfile ~/midway/genos/",data_source,"/plink/chr",chrom,"_rsids --pheno gene.pheno --make-bed --out ",gene," --keep gene.pheno --chr ",chrom,
           " --from-bp ",start," --to-bp ",stop,sep="")
system(cmd)
cmd<-paste("plink --bfile ",gene," --recode tab --out ",gene,sep='')
system(cmd)


#Count number of subjects per allele
count_genos<-function(gene_dat,rsid) {
  plink_rsids<-fread(paste(gene,".map",sep=''))
  line<-which(plink_rsids$V2==rsid)
  plink_genos<-fread(paste(gene,".ped",sep=''))
  IDs<-plink_genos[,1]
  plink_genos<-plink_genos[,7:ncol(plink_genos)]
  genos<-as.data.frame(plink_genos)[,line]
  genos<-gsub(" ", "/", genos)
  genos<-cbind(IDs,genos)
  allele_counts<-matrix(nrow=2,ncol=0)
  for (c in c("A","C","T","G")) {
    num<-sum(str_count(genos$genos,c))
    allele_counts<-cbind(allele_counts,as.data.frame(c(c,num),stringsAsFactors = FALSE))
  }
  allele_counts<-as.data.frame(t(allele_counts),stringsAsFactors = FALSE) %>% arrange(.,V2)
  write.table(t(allele_counts[3:4,]),file=paste(gene,"_",rsid,"_report",sep=''),quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE)
  geno_exp<-merge(gene_dat,genos,by.x="P1",by.y="V1")
  geno_exp$Expression<-as.numeric(as.character(geno_exp$Expression))
  return(geno_exp)
}
geno_exp<-count_genos(gene_dat,rsid)
geno_exp<-geno_exp[geno_exp$genos!="0/0",]

##Add sex
females<-fread(paste("~/midway/expression/",trait,"_F/",tissue,"/",tissue,"_subject_IDs",sep=''),header=FALSE)
females=females$V1
males<-fread(paste("~/midway/expression/",trait,"_M/",tissue,"/",tissue,"_subject_IDs",sep=''),header=FALSE)
males=males$V1

geno_exp$sex=rep(NA,nrow(geno_exp))
geno_exp$sex[geno_exp$P1 %in% females]="female"
geno_exp$sex[geno_exp$P1 %in% males]="male"
geno_exp<-geno_exp[!(is.na(geno_exp$sex)),]

## Boxplot with sex as grouped var
CC_p<-t.test(geno_exp$Expression[geno_exp$genos=="C/C" & geno_exp$sex=="female"],geno_exp$Expression[geno_exp$genos=="C/C" & geno_exp$sex=="male"])$p.value
CA_p<-t.test(geno_exp$Expression[geno_exp$genos=="C/A" & geno_exp$sex=="female"],geno_exp$Expression[geno_exp$genos=="C/A" & geno_exp$sex=="male"])$p.value
AA_p<-t.test(geno_exp$Expression[geno_exp$genos=="A/A" & geno_exp$sex=="female"],geno_exp$Expression[geno_exp$genos=="A/A" & geno_exp$sex=="male"])$p.value
CC_FM_FC=mean(geno_exp$Expression[geno_exp$genos=="C/C" & geno_exp$sex=="female"])/mean(geno_exp$Expression[geno_exp$genos=="C/C" & geno_exp$sex=="male"])
CA_FM_FC=mean(geno_exp$Expression[geno_exp$genos=="C/A" & geno_exp$sex=="female"])/mean(geno_exp$Expression[geno_exp$genos=="C/A" & geno_exp$sex=="male"])
AA_FM_FC=mean(geno_exp$Expression[geno_exp$genos=="A/A" & geno_exp$sex=="female"])/mean(geno_exp$Expression[geno_exp$genos=="A/A" & geno_exp$sex=="male"])

ages<-fread("~/midway/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt")
geno_exp<-merge(ages,geno_exp,by.x="SUBJID",by.y="P1")

B<-ggplot(geno_exp,aes(x=genos,y=Expression,color=sex,fill=genos))+
  geom_boxplot(width=0.5,outlier.shape = NA)+
  geom_jitter(shape=16,size=1.5,position=position_jitterdodge(dodge.width=0.5,jitter.width=0.1),aes(group=sex))+
  ylab("Expression (tpm)")+
  xlab(paste(rsid," genotype",sep=''))+
  scale_y_continuous(limits=c(0,max(geno_exp$Expression+5)))+
  scale_color_manual(values=c("black","black"))+
  scale_fill_manual(values=c(rgb(pony_colors[8,1:3]),rgb(pony_colors[10,1:3]),rgb(pony_colors[13,1:3])))+
  guides(fill=FALSE,color=FALSE)+
  theme_minimal()+
  theme(axis.title=element_text(size=16),
        axis.text=element_text(size=16))

pdf(paste(gene,"_",rsid,"_groupedsex_expression_boxplot.pdf",sep=''),width=8,height=4)
B
dev.off()


##Female boxplot:
F_geno_exp<-geno_exp[geno_exp$sex=="female",] 
F<-ggplot(F_geno_exp,aes(x=genos,y=Expression,fill=genos))+
  geom_boxplot(width=0.5,outlier.shape = NA)+
  geom_jitter(shape=16,size=1.5,position=position_jitter(0.1))+
  ylab("Expression (tpm)")+
  xlab(paste(rsid," genotype",sep=''))+
  scale_y_continuous(limits=c(0,max(geno_exp$Expression+5)))+
  scale_fill_manual(values=c(rgb(pony_colors[8,1:3]),rgb(pony_colors[10,1:3]),rgb(pony_colors[13,1:3])))+
  guides(fill=FALSE)+
  theme_minimal()+
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=18))
F

##Male boxplot:
M_geno_exp<-geno_exp[geno_exp$sex=="male",] 
M<-ggplot(M_geno_exp,aes(x=genos,y=Expression,fill=genos))+
  geom_boxplot(width=0.5,outlier.shape = NA)+
  geom_jitter(shape=16,size=1.5,position=position_jitter(0.1))+
  ylab("Expression (tpm)")+
  xlab(paste(rsid," genotype",sep=''))+
  scale_y_continuous(limits=c(0,max(geno_exp$Expression+5)))+
  scale_fill_manual(values=c(rgb(pony_colors[8,1:3]),rgb(pony_colors[10,1:3]),rgb(pony_colors[13,1:3])))+
  guides(fill=FALSE)+
  theme_minimal()+
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=18))

pdf(paste(gene,"_",rsid,"_sex_expression_boxplot.pdf",sep=''),width=8,height=4)
grid.arrange(F,M,ncol=2)
dev.off()

system(paste("rm gene*pheno ",gene,".*",sep=''))

