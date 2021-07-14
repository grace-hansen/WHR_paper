library(biomaRt)
library(Gviz)

bm <- useMart(host="ensembl.org", biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")

#Tracks
biomTrack <- BiomartGeneRegionTrack(genome="hg38", chromosome=7,start=26290000, end=26420000,name="SNX10 locus", biomart=bm,filter=list(with_refseq_mrna=TRUE))
gtrack <- GenomeAxisTrack(cex=1.5)

pdf("~/papers/TWAS/adipocyte_profiler/ABC/SNX10_locus_genes.pdf",width=7,height=1.5)
plotTracks(list(biomTrack,gtrack),collapseTranscripts="meta",transcriptAnnotation="symbol",
           stackHeight=0.3,
           from=26270000,to=26420000,
           fontsize.group=16,fontcolor.group="black",col="black",fill="black")
dev.off()

