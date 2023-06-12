# Script to get the input parameters and split the gene list for parallelization

# load libraries (these can be installed in the Dockerfile)
library("class")
library("cluster")
library("rpart")
library("lattice")

args <- commandArgs(trailingOnly = TRUE)

# find scripts to source
source_folder <- args[1]

# set directed boolean for downstream analysis
if (args[2] == "Directed") {
  directed <- TRUE
} else {
  directed <- FALSE
}
networkFile = args[3] #ParentChild, BN_digraph_pruned_formatted, (V1(parent), V2(child), tab delim)
targetFile = args[4] #Either local or global seeding gene list. (V1(gene name of network space), V2(group name string)) Break this into parallel jobs. If global, just use all genes in the network in V1, call V2 "Allnodes" or "global". 
outputDirectory = args[5] #where you want this to be saved. This creates a subfolder called KDA. 
number_of_layers = args[6] #normally 6 steps away. 

##Location of KDA R functions, make sure / at the end

# get variables ready for script input
fcausalnet <- networkFile
finputlist <- targetFile
layer <- as.numeric(number_of_layers)
outputDir <- outputDirectory
fgeneinfo <- NULL

# create an output directory
if (directed) {
  dir.create(paste(outputDir, "KeyDriversDirected/", sep = ""))
} else {
  dir.create( paste(outputDir,"KeyDriversUndirected/",sep=""))
}

# source all necessary files
for (f in list.files(source_folder,pattern="*.R")) {
  try(source(paste(source_folder,f,sep=""),local=FALSE))
}

# read in network
cnet <- read.delim(fcausalnet, sep = "\t", header = F)
cnet <- as.matrix(cnet)
cnet = cnet[cnet[,2] != "", ]
totalnodes <- union( cnet[,1] , cnet[,2] )

# get the filename
fname <- gsub(".*/","",getFileName(fcausalnet))

if ( directed ){
  fname <- paste( outputDir,"KeyDriversDirected/", fname , "_L" , layer , sep = "" )
}else{
  fname <- paste( outputDir,"KeyDriversUndirected/", fname, "_L", layer , sep = "" )
}

# read in gene lists
listMatrix <- read.delim(finputlist, sep="\t", header = TRUE)
# dim( listMatrix )
listMatrix <- as.matrix(listMatrix)
# listMatrix[1:2,]
ncols <- dim(listMatrix)[2]
modules <- names(table(listMatrix[,ncols]))
xkdFall <- paste(fname, "_KDx_combined", format(Sys.time(), "%Y-%m-%d %I-%p"), ".txt", sep = "")
xkdFpara <- paste(fname, "_KDx_parameters", format(Sys.time(), "%Y-%m-%d %I-%p"), ".txt", sep = "")
xkdrMatrix <- NULL
paraMatrix <- NULL