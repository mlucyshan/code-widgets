suppressPackageStartupMessages(library("DEXSeq"))
suppressPackageStartupMessages(library("matrixStats"))

inDir=argv[1]
inDir=getwd()
countFiles_exon=list.files(inDir, 
                      pattern="*counts.txt$", 
                      full.names=TRUE)
flattenedFile=argv[2]
sampleTable_exon=data.frame(
  row.names=c("min00_rep1", "min00_rep2", 
              "min30_rep1", "min30_rep2", 
              "min60_rep1", "min60_rep2"),
  condition=c("min00", "min00", "min30", "min30", "min60", "min60"),
  libType=rep("paired-end", 6)
)

sampleTable_exon_subset=data.frame(
  row.names=c("min00_rep1", "min00_rep2", 
              "min30_rep1", "min30_rep2"),
  condition=c("min00", "min00", "min30", "min30"),
  libType=rep("paired-end", 4)
)

dxd_exon_subset=DEXSeqDataSetFromHTSeq(
  countFiles_exon_subset, 
  sampleData=sampleTable_exon_subset, 
  design= ~ sample + exon + condition:exon, 
  flattenedfile=flattenedFile
)
dxd_exon_subset=estimateSizeFactors(dxd_exon_subset)
dxd_exon_subset=estimateDispersions(dxd_exon_subset)
dxd_exon_subset=testForDEU(dxd_exon_subset)
dxd_exon_subset=estimateExonFoldChanges(dxd_exon_subset, fitExpToVar="condition")
dxr_exon_subset=DEXSeqResults(dxd_exon_subset)
table(dxd_exon_subset$pvalue<0.05)
true_dxr_exon_subset=subset(dxd_exon_subset,dxd_exon_subset$pvalue<0.05)
true_groupID_exon_subset=true_dxr2_exon_subset$groupID
pdf(file="/Data04/mengge/moss_mRNA_rai1/exon_usage.pdf", width=11,height=8)
for (i in 1:length(true_groupID_exon)){
  plotDEXSeq(dxr_exon_subset, true_groupID_exon[i], displayTranscripts=TRUE, legend=TRUE,cex.axis=1.2, cex=1.3, lwd=2 )
  plotDEXSeq(dxr_exon_subset, true_groupID_exon[i], legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
}
dev.off()
