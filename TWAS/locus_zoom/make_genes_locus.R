library(biomaRt)
library(Gviz)


bm <- useMart(host="ensembl.org", biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")

#Tracks
biomTrack <- BiomartGeneRegionTrack(genome="hg38", chromosome=7,start=26159590, end=26557590,name="SNX10 locus", biomart=bm,filter=list(with_refseq_mrna=TRUE))
gtrack <- GenomeAxisTrack(cex=1.5)

pdf("~/papers/TWAS/sex_dimorphism/locus_zoom/SNX10_locus_genes.pdf",width=10,height=1.5)
plotTracks(list(biomTrack,gtrack),collapseTranscripts="meta",transcriptAnnotation="symbol",
           stackHeight=0.3,
           start=26159590, end=26557590,
           fontsize.group=18,fontcolor.group="black",col="black",fill="black")
dev.off()

