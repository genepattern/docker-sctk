## The Regents of the University of California and The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2018) by the
## Regents of the University of California abd the 
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

# Load any packages used to in our code to interface with GenePattern.
# Note the use of suppressMessages and suppressWarnings here.  The package
# loading process is often noisy on stderr, which will (by default) cause
# GenePattern to flag the job as failing even when nothing went wrong. 
suppressMessages(suppressWarnings(library(getopt)))
suppressMessages(suppressWarnings(library(optparse)))

suppressPackageStartupMessages({
  library(MAST)
  library(singleCellTK)
  library(xtable)
})


# Print the sessionInfo so that there is a listing of loaded packages, 
# the current version of R, and other environmental information in our
# stdout file.  This can be useful for reproducibility, troubleshooting
# and comparing between runs.
sessionInfo()

# Get the command line arguments.  We'll process these with optparse.
# https://cran.r-project.org/web/packages/optparse/index.html
arguments <- commandArgs(trailingOnly=TRUE)

print(packageVersion("singleCellTK"))
# Declare an option list for optparse to use in parsing the command line.
option_list <- list(
  # Note: it's not necessary for the names to match here, it's just a convention
  # to keep things consistent.
  make_option("--assay.file", dest="assay.file"),
  make_option("--assay.name", dest="assay.name"),
  make_option("--cls.file",  dest="cls.file"),
  make_option("--Run.PCA",  dest="Run.PCA", type="logical"),
  make_option("--Run.TSNE", dest="Run.TSNE", type="logical"),
  make_option("--output.file", dest="output.file"),
  make_option("--log.transform", dest="log.transform", type="logical")

  )

# Parse the command line arguments with the option list, printing the result
# to give a record as with sessionInfo.
opt <- parse_args(OptionParser(option_list=option_list),  args=arguments)
print("Options are ")
print(opt)

# Load some common GP utility code for handling GCT files and so on.  This is included
# with the module and so it will be found in the same location as this script (libdir).
source(file.path("/usr/local/bin/sctk/", "common.R"))

# Process the parameters.  
# Note that since some parameters are optional, we must be prepared to receive no value
# Also check for blank values since these are ignored.
if (is.null(opt$assay.file) || grepl("^[[:space:]]*$", opt$assay.file)) {
   assay.file <- NULL
} else {
   assay.file <- opt$assay.file
}

if (is.null(opt$assay.name) || grepl("^[[:space:]]*$", opt$assay.name)) {
   assay.name <- "assay"
} else {
   assay.name <- opt$assay.name
}

# read the data frame and get the data
gct <-read.gct(assay.file)
df = data.frame(gct[2])

if (opt$log.transform) {
	df = log2(df)
}

if (is.null(opt$cls.file) || grepl("^[[:space:]]*$", opt$cls.file)) {
   stop("Required parameter cls file was not provided.")
} else {
   cls <- read.cls(opt$cls.file)
   cls.label <- c("sample","class")
   cdf = data.frame(colnames(df), cls$labels)
   rownames(cdf) <- c(colnames(df))
   condition="cls.labels" 
}


gct_sce <- createSCE(assayFile=df, annotFile=cdf, assayName=assay.name, inputDataFrames=TRUE)

if (opt$Run.PCA){
	gct_sce <- getPCA(gct_sce, useAssay=assay.name, reducedDimName="PCA")
	pdf(paste(opt$output.file, "_PCA.pdf", sep=""))
	print(plotPCA(gct_sce, reducedDimName = "PCA", colorBy = condition))
	dev.off()
	gctPCA = {}
	gctPCA$data <-data.frame(reducedDim(gct_sce, "PCA"))
	write.gct(gctPCA, file.path(getwd(), paste(opt$output.file, "_PCA.gct", sep="")))
}

if (opt$Run.TSNE){
	gct_sce <- getTSNE(gct_sce, useAssay=assay.name, reducedDimName="TSNE")
	pdf(paste(opt$output.file, "_TSNE.pdf", sep=""))
	print(plotTSNE(gct_sce, reducedDimName = "TSNE", colorBy = condition))
	dev.off()
	gctTSNE = {}
	gctTSNE$data <-data.frame(reducedDim(gct_sce, "TSNE"))
	write.gct(gctTSNE, file.path(getwd(), paste(opt$output.file, "_TSNE.gct", sep="")))
}

mast_results <- MAST(gct_sce, condition = condition, useThresh = TRUE, useAssay = assay.name)
pdf(paste(opt$output.file, "_Violin.pdf", sep=""))
print(MASTviolin(gct_sce, useAssay = assay.name , fcHurdleSig = mast_results, threshP = TRUE, condition = condition))
dev.off()

pdf(paste(opt$output.file, "_Regression.pdf", sep=""))
print(MASTregression(gct_sce, useAssay = assay.name, fcHurdleSig = mast_results, threshP = TRUE, condition = condition))
dev.off()

pdf(paste(opt$output.file, "_DiffEx.pdf", sep=""))
print(plotDiffEx(gct_sce, useAssay = assay.name, condition = condition, geneList = mast_results$Gene[1:100], annotationColors = "auto",  clusterRow=FALSE, displayRowLabels = FALSE, displayColumnLabels = FALSE))
dev.off()






