## try http:// if https:// URLs are not supported
# module specific packages first 

source("http://bioconductor.org/biocLite.R")
biocLite("MAST")
biocLite("singleCellTK")
biocLite("xtable")
biocLite("bladderbatch")

suppressPackageStartupMessages({
  library(MAST)
  library(singleCellTK)
  library(xtable)
})

