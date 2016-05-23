#!/bin/Rscript 
args<-commandArgs(TRUE)
inputDir<-args[1]
outputDirMaster<-args[2]
#check if SeqLogo is available
if ("seqLogo" %in% rownames(installed.packages())==FALSE){
  source("http://bioconductor.org/biocLite.R")
  biocLite("seqLogo")
}
#import the SeqLogo library
require("seqLogo")
require("tools")

#Create the directory structure
inputDir<-paste(inputDir, "homerResults/", sep="/")
print(inputDir)
#inputDir<-"/Users/mshan/Desktop/HOMER_output_moss_xlink_00min_CSAR_PPS_FDR05/homerResults/"
expName<-strsplit(inputDir, split="/")
expName<-expName[[1]][length(expName[[1]])-1]
#outputDirMaster<-"/Users/mshan/Desktop/HOMER_seqLogos2"
outputDir<-paste(outputDirMaster, expName, sep="/")
dir.create(outputDir, showWarnings=TRUE, recursive=TRUE, mode="0777")
#Read in the HOMER output motif files
files<-dir(inputDir, pattern="motif[0-9]+\\.motif")
for (f in files){
  motifName<-tools::file_path_sans_ext(basename(f))
  titleLine<-strsplit((readLines(paste(inputDir, f, sep="/"), n=1)), split="\t")
  titleLine<-strsplit((titleLine[[1]][length(titleLine[[1]])]), split=",")
  pvalue<-titleLine[[1]][length(titleLine[[1]])]
  weightMatrix<-t(read.table(paste(inputDir, f, sep="/"), skip=1))
  pwm<-makePWM(weightMatrix)
  outputPlotName<-paste(outputDir, paste(paste(motifName, "seqLogo", sep="_"),"png", sep="."), sep="/")
  print(outputPlotName)
  png(outputPlotName)
  par(oma=c(4,1,1,1))
  plot.new()
  seqLogo(pwm)
  title(motifName)
  par(fig=c(0,1,0,1), oma=c(0,0,0,0), mar=c(0,0,0,0), new=TRUE)
  plot(0,0, type="n", bty="n", xaxt="n", yaxt="n")
  legend(x="topleft", legend=c(pvalue), cex=1.5, bty="n", xpd=TRUE, inset=c(0,0))
  dev.off()
}
