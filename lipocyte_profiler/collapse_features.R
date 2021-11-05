#####remove redundant features #####
####keep redundant feature with largest effect size#####
####input data you want to filter####
####eg. df = significant features##
library(data.table)
library(tidyverse)
library(stringr)

collapse_features<-function(df) {
  all_features<-fread(df) #All features, as defined by Sam and Sophie
  all_features$collapsed_features<-all_features$features
  
  ######## Give redundant features same name #######
  all_features$collapsed_features<-str_remove(all_features$collapsed_features, "[0-9][0-9][_][0-9][0-9]")
  all_features$collapsed_features<-str_remove(all_features$collapsed_features, "[0-9][_][0-9][0-9]")
  all_features$collapsed_features<-str_remove(all_features$collapsed_features, "[0-9][_][0-9]")
  all_features$collapsed_features<-str_remove(all_features$collapsed_features, "[_][0-9][_]")
  all_features$collapsed_features<-str_remove(all_features$collapsed_features, "[_][0-9][0-9][_]")
  all_features$collapsed_features<-str_remove(all_features$collapsed_features, "[0-9]*[0-9]")
  
  ####### Order features based on feature string and effect size ######
  all_features$SNP.eta_sq<-as.numeric(all_features$SNP.eta_sq)
  all_features<-all_features[order(all_features$collapsed_features,-abs(all_features$SNP.eta_sq)),]
  
  ####### Keep redundant feature with highest effect size
  nonredundant_all_features<-all_features[!duplicated(all_features$collapsed_features), ]
  return(nonredundant_all_features)
}

for (tp in c("day0","day3","day8","day14")) {
  for (ct in c("sc","vc")) {
    df<-collapse_features(paste("~/medusa/papers/TWAS/lipocyte_profiler/scatter/rs1534696_allfeatures_",tp,"_",ct,"_bothsexes.tsv",sep=''))
    write.table(df,paste("~/medusa/papers/TWAS/lipocyte_profiler/scatter/rs1534696_allfeatures_",tp,"_",ct,"_bothsexes_nonredundant.txt",sep=''),sep='\t',quote=FALSE,row.names=FALSE)
  }
}

######### Evaluate whether the same features are sig in this collapse and the original collapse #####
#sig_features<-fread("~/medusa/papers/TWAS/lipocyte_profiler/scatter/rs1534696_sigfeatures_day14_bothsexes.tsv") #Significant features, as defined by Sam and Sophie
#sig_features$collapsed_features<-sig_features$raw_features
#sig_features$collapsed_features<-str_remove(sig_features$collapsed_features, "[0-9][0-9][_][0-9][0-9]")
#sig_features$collapsed_features<-str_remove(sig_features$collapsed_features, "[0-9][_][0-9][0-9]")
#sig_features$collapsed_features<-str_remove(sig_features$collapsed_features, "[0-9][_][0-9]")
#sig_features$collapsed_features<-str_remove(sig_features$collapsed_features, "[_][0-9][_]")
#sig_features$collapsed_features<-str_remove(sig_features$collapsed_features, "[0-9]*[0-9]")

#nonredundant_all_features<-collapse_features("~/medusa/papers/TWAS/lipocyte_profiler/scatter/rs1534696_allfeatures_day14_sc_bothsexes.tsv")
#sig_features_GH<-nonredundant_all_features[nonredundant_all_features$q<0.1]
#sum(!(sig_features$collapsed_features %in% sig_features_GH$collapsed_features)) #Should equal 0: all sig features should be in both lists
