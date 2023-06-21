#!/usr/bin/env Rscript
#
##-----------------------------------------------------------------------
# Input:
#   1. directed or undirected network: "Directed" OR "Undirected";
#   2. causal network: #ParentChild, BN_digraph_pruned_formatted,
# (V1(parent), V2(child), tab delim);
#   3. gene list: #Either local or global seeding gene list.
# (V1(gene name of network space), V2(group name string))
#   Break this into parallel jobs.
#   If global, just use all genes in the network in V1,
#   call V2 "Allnodes" or "global".
#   4. outputDirectory: #where you want this to be saved.
#   5. expand subnetwork based on L-layer neighbors: layer=0
#   NOTUSED. gene annotation file(NULL if not available): fgeneinfo
#
# Output:
#    1. keydrivers for each subnetwork: "*_keydriver.txt"
#    2. Cytoscape network: *_cys.txt #Troubleshooting
#    3. Cytoscape node properties: *_cys-nodes.txt #Troubleshooting
#    4. combined keydrivers for all runs: "_KDx_combined.txt"
# Likely turning this off to make your own combined
#    5. parameters used for key driver analysis: "_KDx_parameters.txt"
#
#
# -------------------- Parameters to be changed -------------------------
library("class")
library("cluster")
library("rpart")
library("lattice") # require is design for use inside functions 
library("glue")

args <- commandArgs(trailingOnly = TRUE)

if (args[1] %in% c("directed")) {
  directed <- TRUE
} else {
  directed <- FALSE
}
networkFile = args[2] #ParentChild, BN_digraph_pruned_formatted, (V1(parent), V2(child), tab delim)
targetFile = args[3] #Either local or global seeding gene list. (V1(gene name of network space), V2(group name string)) Break this into parallel jobs. If global, just use all genes in the network in V1, call V2 "Allnodes" or "global". 
#outputDirectory = args[5] #where you want this to be saved. This creates a subfolder called KDA. 
number_of_layers = args[4] #normally 6 steps away. 

functions <- args[5] 
source(glue("./{functions}")

##Location of KDA R functions, make sure / at the end
# source functions:
#source(glue("{source_folder}/concatenate.R"))
#source(glue("{source_folder}/configureNodeVisualization.R"))
#source(glue("{source_folder}/degreeByLinkPairs.R"))
#source(glue("{source_folder}/downStreamGenes.R"))
#source(glue("{source_folder}/findNLayerNeighborsLinkPairs.R"))
#source(glue("{source_folder}/getAllParts.R"))
#source(glue("{source_folder}/getFileName.R"))
#source(glue("{source_folder}/getFileExtension.R"))
#source(glue("{source_folder}/getMatchedIndexFast.R"))
#source(glue("{source_folder}/getSubnetworkLinkPairs.R"))
#source(glue("{source_folder}/keyDriverAnalysis.R"))
#source(glue("{source_folder}/keydriverInSubnetwork.R"))
#source(glue("{source_folder}/makeSNP.R"))
#source(glue("{source_folder}/mergeTwoMatricesByKeepAllPrimary.R"))
#source(glue("{source_folder}/removeDuplicatedLinks.R"))
#source(glue("{source_folder}/replaceString.R"))
#source(glue("{source_folder}/setElementInSet.R"))
#source(glue("{source_folder}/setInSets.R"))
#source(glue("{source_folder}/setsub.R"))
#source(glue("{source_folder}/splitString.R"))

#KDARfunctions = source_folder

fcausalnet <- networkFile
finputlist <- targetFile
layer <- as.numeric(number_of_layers)
#outputDir <- outputDirectory
fgeneinfo <- NULL

# 2. specify the directory for holding analysis results

#if (directed) {
#  dir.create(paste(outputDir, "KeyDriversDirected/", sep = ""))
#} else {
#  dir.create( paste(outputDir,"KeyDriversUndirected/", sep=""))
#}

# -----------------------------End of Parameters to be changed --------------------------------------
#install.packages("class")
#install.packages("cluster")
#install.packages("rpart")
#install.packages("lattice")
#install.packages("SpaDES")


# Windows
# memory.size( TRUE )   # check the maximum memory that can be allocated
# memory.limit( size = 3800 )   # increase the available memory

# box::use(workspace/kda_nextflow/lib/getFileName.R[...])

# source scripts with path to script folder
#for (f in list.files(KDARfunctions,pattern="*.R")) {
#  try(source(glue(KDARfunctions,f)))
#}


################################################################################################
#    1. read in network

cnet <- read.delim(fcausalnet, sep = "\t", header = F)
cnet <- as.matrix(cnet)
cnet = cnet[cnet[,2] != "", ]

totalnodes <- union( cnet[,1] , cnet[,2] )

fname <- gsub(".*/","",getFileName(fcausalnet))
#fname <- gsub(".*/","",fcausalnet)
print(fname)

#if ( directed ){
#  fname <- paste( outputDir,"KeyDriversDirected/", fname , "_L" , layer , sep = "" )
#}else{
#  fname <- paste( outputDir,"KeyDriversUndirected/", fname, "_L", layer , 
#  sep = "" )
#}

fname <- paste(fname , "_L" , layer , sep = "")

################################################################################
# 2. read in gene lists

listMatrix <- read.delim( finputlist , sep="\t" , header = TRUE )
# dim( listMatrix )
listMatrix <- as.matrix( listMatrix )
# listMatrix[1:2,]
ncols <- dim( listMatrix )[2]

modules <- names( table( listMatrix[,ncols] ) )

#xkdFall <- paste( fname , "_KDx_combined.txt" , sep = "" )
#xkdFpara <- paste( fname , "_KDx_parameters.txt" , sep = "" )
xkdFall <- paste( fname , "_KDx_combined", format(Sys.time(), "%Y-%m-%d %I-%p"), ".txt" , sep = "" )
xkdFpara <- paste( fname , "_KDx_parameters", format(Sys.time(), "%Y-%m-%d %I-%p"), ".txt" , sep = "" )
xkdrMatrix <- NULL
paraMatrix <- NULL

############################################################################ 
# 3. process each gene list
#
#em=modules[1]
for ( em in modules ){
  print( paste( "*****************" , em , "********************" ) )
  esel <- listMatrix[,ncols] == em
  # remove abnormal gene names
  genes <- union( listMatrix[esel,1] , NULL )
  genes <- genes[genes != ""]
  genes <- genes[!is.na( genes )]
  no.genes <- length( genes )
  
  em2 <- gsub( em , ":" , "" )
  em2 <- gsub( em2 , " " , "-" )
  
  key2 <- paste( fname , "_KD_" , em2 , sep = "" )
  onetFname <- paste( key2 , ".pair" , sep = "" )
  snpFname <- paste( key2 , ".snp" , sep = "" )
  kdFname <- paste( key2 , "_keydriver.txt" , sep = "" )
  
  if(layer >=1 ){
     # expand network by K-hop nearest neighbors layers
     expandNet <- findNLayerNeighborsLinkPairs( linkpairs = cnet , subnetNodes = genes ,
			                                     nlayers = layer , directed = FALSE )
  }else{
   # no expansion
     expandNet <- getSubnetworkLinkPairs( linkpairs = cnet , subnetNodes = genes )
  }
  dim( expandNet )
  
  allnodes <- union( expandNet[,1] , expandNet[,2] )
  #write.table(expandNet, onetFname, sep="\t",quote=FALSE, col.names=T, row.names=FALSE)

  ################################################################################################
  # 4. keydriver for a given network
  #
  if (directed){
    ret <- keydriverInSubnetwork( linkpairs = expandNet , signature = genes, background=NULL, directed = directed ,
			         nlayers = layer , enrichedNodes_percent_cut=-1, FET_pvalue_cut=0.05,
			         boost_hubs=T, dynamic_search=T, bonferroni_correction=T, expanded_network_as_signature =F)
  }else{
    ret <- keydriverInSubnetwork( linkpairs = expandNet , signature = genes , directed = directed ,
			         nlayers = layer , enrichedNodes_percent_cut=-1, FET_pvalue_cut=0.05,
				boost_hubs=T, dynamic_search=T, bonferroni_correction=T, expanded_network_as_signature =F)
  }
  if ( is.null( ret ) ){
  	next
  }
  
  fkd <- ret[[1]]
  parameters <- ret[[2]]
  
  fkd2 <- cbind( rep( em , dim( fkd )[1] ) , fkd )
  xkdrMatrix <- rbind( xkdrMatrix , fkd2 )
  
  paraMatrix <- rbind( paraMatrix , c( key2 , parameters ) )
  
  write.table( fkd , kdFname , sep = "\t" , quote = FALSE , col.names = TRUE , row.names = FALSE )

}


# save all key drivers
colnames( xkdrMatrix ) <- c( "module" , colnames( fkd ) )
write.table( xkdrMatrix , xkdFall , sep = "\t" , quote = FALSE , col.names = TRUE , row.names = FALSE )

# save parameters used
#
colnames( paraMatrix ) <- c( "subnet" , colnames( parameters ) )
write.table( paraMatrix , xkdFpara , sep = "\t" , quote = FALSE , col.names = TRUE ,
		     row.names = FALSE )

if( !is.null( fgeneinfo ) ){
   infoMatrix <- read.delim( fgeneinfo , sep = "\t" , header = TRUE )
   dim( infoMatrix )
   xkdrMatrix2 <- cbind( xkdrMatrix , c( 1:( dim( xkdrMatrix )[1] ) ) )

   ic <- dim( infoMatrix )[2] + 1
   merged <- merge( infoMatrix , xkdrMatrix2 , by.x = 1 , by.y = 2 , all.y = TRUE )
   merged <- as.matrix( merged )
   ic2 <- dim( merged )[2]
   xf <- merged[,c( ic , setdiff( 1:ic2 , ic ) )]
   ic3 <- dim( merged )[2]
   mo <- order( as.integer( merged[,ic3] ) )
   write.table( xf[mo,-ic3] , xkdFall , sep = "\t" , quote = FALSE , col.names = TRUE , row.names = FALSE )
}

## ------------------------------------- END ------------------------------------------------
