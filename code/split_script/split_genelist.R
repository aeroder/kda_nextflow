# Script to get the input parameters and split the gene list for parallelization

# load libraries (these can be installed in the Dockerfile)
library("class")
library("cluster")
library("rpart")
library("lattice")

# gather the arguments, should only be one
args <- commandArgs(trailingOnly = TRUE)
if(!length(args) == 1) {
  stop("gene list only should be provided as two column,
  tab separated file")
}

# get the gene list name
finputlist = args[1]  

finputlist <- "~/workspace/kda_nextflow/input/KDAinputfile.txt"


# read in gene lists
listMatrix <- read.delim(finputlist, sep="\t", header = FALSE) |> as.matrix()
rownames(listMatrix) <- NULL
# listMatrix[1:2,]
ncols <- dim(listMatrix)[2]
modules <- names(table(listMatrix[,ncols]))