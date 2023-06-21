#!/usr/bin/env Rscript
concatenate <- function( myvect , mysep="" ){
  if ( is.null( myvect ) )
  {
    return ( "" )
  }
  else if ( isTRUE( all.equal( length( myvect ) , 1 ) ) )
  {
    return ( as.character( myvect ) )
  }

  return( paste( as.character( myvect ), sep = "" , collapse = mysep ) )

}

configureNodeVisualization = function(allnodes, signature, kdaMatrix, bNodeSz=40, bFontSz=12) {

   # SIG--signature; NSIG--not signature; GKD--Global KeyDriver; LKD--Local KeyDriver; NKD--Not KeyDriver
   #
   xcategories = c("SIG_GKD", "SIG_LKD", "SIG_NKD", "NSIG_GKD", "NSIG_LKD", "NSIG_NKD"); xcat2=c("NSIG_GKD", "NSIG_LKD", "NSIG_NKD")
   xcolors     = c("red",     "blue",    "lightgreen",  "red",      "blue",     "grey");  names(xcolors)<- xcategories
   xshapes     = c("square",  "square",  "circle",   "circle",    "circle",   "circle");names(xshapes)<- xcategories
   xsizes      = c(3*bNodeSz, 2*bNodeSz, bNodeSz,   3*bNodeSz,  2*bNodeSz,  bNodeSz);     names(xsizes) <- xcategories
   xfontsz     = c(3*bFontSz, 2*bFontSz, bFontSz,   3*bFontSz,  2*bFontSz,  bFontSz);     names(xfontsz)<- xcategories

   no.nodes = length(allnodes)

   # legend table 
   legendtb = cbind(xcategories, xshapes, xcolors, xcolors, xsizes, xfontsz)
   colnames(legendtb) <- c("label", "shape", "color", "border", "node_size", "font_size")

   sigInNet = intersect(allnodes, signature)
   sig_status = rep("NSIG", no.nodes); names(sig_status) <- allnodes; sig_status[sigInNet]="SIG"
   kdr_status = rep("NKD",  no.nodes); names(kdr_status) <- allnodes; 

   nf.cols = dim(kdaMatrix)[2]; nf.rows = dim(kdaMatrix)[1]
   keydrvNames = NULL
   if(nf.rows>0) {
     keydrv  = as.integer(kdaMatrix[,nf.cols])

     # global driver
     keysel  = c(1:nf.rows)[keydrv==1]; 
     keydrvNames = kdaMatrix[keysel,1];
     kdr_status[keydrvNames] = "GKD"

     # local driver
     if(sum(keydrv==0)>0) {
        keysel  = c(1:nf.rows)[keydrv==0];
        keydrvNames = kdaMatrix[keysel,1];
        kdr_status[keydrvNames] = "LKD"
     }

     # combined signature-keydriver status
     #
     sigkdr_status=paste(sig_status, kdr_status, sep="_")
     hnList = tapply(allnodes, sigkdr_status, list) # make a list for each category
     sigkdr_names = names(hnList)

     isNonSig = intersect(xcat2, sigkdr_names) # if all nodes are signatures, we use only circle for display
     if(length(isNonSig)==0){
       xshapes     = c("circle",   "circle",   "circle",  "circle",   "circle",   "circle");names(xshapes)<- xcategories
     }

     # set up actual visualization properties
     yHighColor = xcolors[sigkdr_names]
     yHighShape = xshapes[sigkdr_names]
     yHighSize  = xsizes[sigkdr_names]
     yHighFontSZ= xfontsz[sigkdr_names]

   } else {
     hnList     = list(sigInNet) # highlight only signature
     yHighColor = c("brown")
     yHighShape = c("circle")
     yHighSize  = c("1")
     yHighFontSZ= c("1")
   }

   return( list(hnList, cbind(yHighColor, yHighShape, yHighSize, yHighFontSZ), legendtb) )
}

degreeByLinkPairs <- function( linkpairs , directed = F , cleangarbage = F )
{
    codepair <- c( 0 , 1 )  #[1] for no connection, [2] for connection

    edgesInNet <- dim( linkpairs )[1]
    
    # consider both columns
	# Need to find out why we build the object this way
	#  Could be a cleaner way, but there may be a reason why
	#  it's currently done this way
    allnodenames <- NULL
    allnodenames <- c( allnodenames , as.character( linkpairs[,1] ) )
    allnodenames <- c( allnodenames , as.character( linkpairs[,2] ) )
     
    nametable <- table( allnodenames )
    # Not sure why this was here, uncommented out
	# length( nametable )
    
    uniquenames <- names( nametable )
    no.uniquenames <- length( uniquenames )

    totallinks <- as.integer( nametable ) # no of links for each node
    totalmatrix <- cbind( names( nametable ),  totallinks )

    if ( directed )
	{
        # outlines
        dnodenames <- as.character( linkpairs[,1] )
        dnametable <- table( dnodenames )
        duniquenames <- names( dnametable )
        dmatrix <- cbind( names( dnametable ), as.integer( dnametable ) )
        colnames( dmatrix ) <- c( "node" , "links" )

        iolinks <- mergeTwoMatricesByKeepAllPrimary( primaryMatrix = cbind( uniquenames ) ,
				minorMatrix = dmatrix , missinglabel = "0" , keepAllPrimary = TRUE ,
				keepPrimaryOrder = TRUE , keepAll = FALSE )
        outlinks <- as.integer( as.matrix( iolinks[,2] ) )


        # inlines
        dnodenames <- as.character( linkpairs[,2] )
        dnametable <- table( dnodenames )
        duniquenames <- names( dnametable )
        dmatrix <- cbind( names( dnametable ) , as.integer( dnametable ) )
        colnames( dmatrix ) <- c( "node" , "links" )

        iolinks <- mergeTwoMatricesByKeepAllPrimary( primaryMatrix = cbind( uniquenames ) ,
				minorMatrix = dmatrix , missinglabel = "0" , keepAllPrimary = TRUE ,
				keepPrimaryOrder = TRUE , keepAll = FALSE )
        inlinks <- as.integer( as.matrix( iolinks[,2] ) ) 

    }
	else
	{
        inlinks <- totallinks
        outlinks <- totallinks
    }

    #hubidx    = order(-totallinks)

    # output in/out links for each gene
    #
    linksMatrix <- cbind( inlinks , outlinks , totallinks )
    colnames( linksMatrix ) <- c( "inlinks" , "outlinks" , "totallinks" )
    rownames( linksMatrix ) <- uniquenames

    rm( inlinks , outlinks , totallinks )

    if ( cleangarbage )
	{
        collect_garbage()
    }

    return( data.frame( linksMatrix ) )
}

downStreamGenes <- function( netpairs , seednodes , N = 100 , directed = TRUE )
{
   prenodes = seednodes
   cnt = N
   while(T) {
      retlinks = findNLayerNeighborsLinkPairs(linkpairs=netpairs, subnetNodes=prenodes, 
                      nlayers=1, directed=directed)

      if(is.null(retlinks)){return (NULL); }

      curnodes = union(retlinks[,1],retlinks[,2]) 
      pcdiff   = setdiff(curnodes, prenodes)
      prenodes = curnodes
      
      if(length(pcdiff)==0){break}

      cnt= cnt-1
      if (cnt==0) {break;}
   }

   if(is.null(retlinks)){
        return (NULL)

   } else {
        return(curnodes)
   }
}

findNLayerNeighborsLinkPairs <- function( linkpairs , subnetNodes , nlayers = 1 , directed = FALSE )
{
  #linkpairs=linkpairs; subnetNodes=ineighbors; nlayers=nlayers-1; directed=directed
   merged <- merge( linkpairs , subnetNodes , by.x = 1 , by.y = 1 , all = FALSE )
   merged <- as.matrix( merged )

   if ( !directed )
   {
# undirected networks
     mergeleft <- merge( linkpairs , subnetNodes , by.x = 2 , by.y = 1 , all = FALSE )
     mergeleft <- as.matrix( mergeleft )
     mergeleft <- mergeleft[,c( 2 , 1 )] # keep the original link direction
     merged <- rbind( merged , mergeleft )
   }

   if ( isTRUE( all.equal( dim( merged )[1] , 0 ) ) )
   {
	   return( NULL )
   }

   if ( isTRUE( all.equal( dim( merged )[1] , 1 ) ) &
	    isTRUE( all.equal( merged[1,1] , merged[1,2] ) ) )
   {
	   return( merged )
   }

   dim1 <- dim( merged )[1]
   if ( isTRUE( all.equal( dim1 , 0 ) ) )
   {
# no links
         return( NULL )
   }
   else if ( is.null( dim1 ) )
   {
# only one link
       merged <- rbind( merged )
   }
   merged <- removeDuplicatedLinks( merged , directed )
   if ( isTRUE( all.equal( nlayers , 1 ) ) )
   {
      return ( merged )
   }

   # get nodes   
   ineighbors <- union( merged[,1] , merged[,2] )

   if (nlayers==1){
      res=getSubnetwork_LinkPairs(linkpairs, subnetNodes=ineighbors)
      return (res)
   }


   # stop earlier if no change
   #
   common <- intersect( ineighbors , subnetNodes )
   if ( length( common ) == length( ineighbors ) )
   {
      return ( merged )
   }

   return ( findNLayerNeighborsLinkPairs( linkpairs , ineighbors , nlayers - 1 , directed ) )
}

getAllParts <- function( fullfnames , sep = "-" , retLen = FALSE )
{
  splitted <- unlist( strsplit( fullfnames[1] , sep ) )
  nn <- length( fullfnames )
  nf <- length( splitted )
  ret <- matrix( "" , nn , nf )
  lens <- rep( 0 , nn )
  for( i in c( 1:nn ) )
  {
    each <- fullfnames[i]
    splitted <- unlist( strsplit( each , sep ) )
    ino <- length( splitted )
    if ( ino >= nf )
	{
       ret[i,] <- splitted[1:nf]
    }
	else
	{
       ret[i,] <- c( splitted , rep( "" , nf - ino ) )
    }
    lens[i] <- ino
  }

  if ( retLen )
  {
     return( lens )
  }
  else
  {
     return( ret ) 
  }
  
}

getFileExtension <- function( fullfname )
{
    splitted <- unlist( strsplit( fullfname , "\\." ) )

    if ( length( splitted ) > 1 )
	{
      return( splitted[length( splitted )] )
    }
	else
	{
      return( "" )
    }
}

getFileName <- function( fullfname )
{
    ext <- getFileExtension( fullfname )
    if (ext == "" )
	{
       return( fullfname )
    }
    extd <- paste( "." , ext , sep = "" )
    return( splitString( fullfname , extd )[1] )
}

getMatchedIndexFast=function(cvector, subvect){
  fullindex = c(1:length(cvector) )
  orgIdx    = cbind(cvector, fullindex)

  index2    = c(1:length(subvect))
  subIdex   = cbind(subvect, index2)

  merged    = merge(subIdex, orgIdx, by.x=1, by.y=1, all.x=T)
  merged    = as.matrix(merged)

  if(dim(merged)[1]>1){
    od        = order(as.integer(merged[,2]))  # restore the original order of subvect
    merged    = merged[od, ]
  }
  
  outIndex  = as.integer(merged[,3])

  return (outIndex)
}

getSubnetworkLinkPairs <- function( linkpairs , subnetNodes )
{
   mergeright <- merge( linkpairs , subnetNodes , by.x = 2 , by.y = 1 , all = FALSE )
   if ( isTRUE( all.equal( dim( mergeright )[1] , 0 ) ) )
   {
     return( NULL )
   }

   mergeleft <- merge( mergeright , subnetNodes , by.x = 2 , by.y = 1 , all = FALSE )

   #mergeright = merge(linkpairs, subnetNodes, by.x=1, by.y=1, all=F)
   #mergeleft2 = merge(mergeright, subnetNodes, by.x=2, by.y=1, all=F)

   if ( isTRUE( all.equal( dim( mergeleft )[1] , 0 ) ) )
   {
     return( NULL )
   }

   return( as.matrix( mergeleft ) )   
}

keyDriverAnalysis = function(inputnetwork, signature, directed=T, nlayer_expansion=1,
      nlayer_search=6, enrichedNodes_percent_cut=-1, boost_hubs=T, dynamic_search=T, 
      FET_pvalue_cut=0.05, use_corrected_pvalue=T, outputfile=NULL, expanded_network_as_signature=FALSE) 
{

   if (!is.null(outputfile)) {
      #onetFname   = paste(outputfile, outputDir, key2, ".pair", sep='')
      snpFname    = paste(outputfile, ".snp",  sep='')
      kdFname     = paste(outputfile, "_keydriver.xls",  sep='')
   }

   # overlap between network & signature
   wholenodes = union(inputnetwork[,1], inputnetwork[,2]); no.wholenodes=length(wholenodes)
   wholeOvlp  = intersect(wholenodes, signature); no.wholeOvlp = length(wholeOvlp)
   
   if(length(wholeOvlp)<=2) {return (NULL)}

   if(nlayer_expansion >=1 ) {
      # expand network by n-layer nearest neighbors
      expandNet = findNLayerNeighborsLinkPairs(linkpairs=inputnetwork, 
                  subnetNodes=signature, nlayers=nlayer_expansion, directed=directed)
   } else if(nlayer_expansion ==0 ){
      # no expansion
      expandNet = getSubnetworkLinkPairs(linkpairs=inputnetwork, subnetNodes=signature)
   } else{
      expandNet = inputnetwork
   }

   if(is.null(expandNet)) {return (NULL)}
   if(dim(expandNet)[1]<=10) {return (NULL)}

   dim(expandNet)
   print(paste("dim(expandNet): ", dim(expandNet)) )

   allnodes = sort(union(expandNet[,1], expandNet[,2])); no.nodes=length(allnodes)

   # convert IDs into indices
   netIdxSrc = getMatchedIndexFast(allnodes, expandNet[,1])
   netIdxDst = getMatchedIndexFast(allnodes, expandNet[,2])
   signatIdx = getMatchedIndexFast(allnodes, intersect(allnodes, signature))
   expandNetIdx = cbind(netIdxSrc, netIdxDst)

   ################################################################################################
   # 4. keydriver for a given network
   #
   #linkpairs=expandNetIdx; signature=signatIdx; background=c(no.wholenodes, no.wholeOvlp); directed=directed; nlayers=nlayer_search; enrichedNodes_percent_cut=enrichedNodes_percent_cut; FET_pvalue_cut=FET_pvalue_cut; boost_hubs=boost_hubs; dynamic_search=dynamic_search; bonferroni_correction=use_corrected_pvalue

   if (directed) {
     ret= keydriverInSubnetwork(linkpairs=expandNetIdx, signature=signatIdx, background=c(no.wholenodes, no.wholeOvlp),
              directed=directed, nlayers=nlayer_search, enrichedNodes_percent_cut=enrichedNodes_percent_cut, 
              FET_pvalue_cut=FET_pvalue_cut, 
              boost_hubs=boost_hubs, dynamic_search=dynamic_search, bonferroni_correction=use_corrected_pvalue, expanded_network_as_signature = expanded_network_as_signature)
   } else{
     ret= keydriverInSubnetwork(linkpairs=expandNetIdx, signature=signatIdx, background=c(no.wholenodes, no.wholeOvlp),
              directed=directed, nlayers=nlayer_search, enrichedNodes_percent_cut=enrichedNodes_percent_cut,
              FET_pvalue_cut=FET_pvalue_cut, 
              boost_hubs=boost_hubs, dynamic_search=dynamic_search, bonferroni_correction=use_corrected_pvalue, expanded_network_as_signature = expanded_network_as_signature)
   }

   if ( is.null(ret)) {return (NULL)}

   # retrieve results
   #
   fkd = ret[[1]]
   parameters = ret[[2]]

   fkd[,1] = allnodes[as.integer(fkd[,1])]   

   nodeDegree = ret[[3]]
   nodeDegree[,1] = allnodes[as.integer(nodeDegree[,1])]   

   if (!is.null(outputfile)) {
      write.table(fkd, kdFname, sep="\t",quote=FALSE, col.names=T, row.names=FALSE)

      ################################################################################################
      #  output networks & key drivers for visualization
      #
      #     Cytoscape output: 1) network file - *_cys.txt 2) node property file: *_cys-nodes.txt
      #  * signature== genes need be corrected in other version
      # allnodes=allnodes; signature=signature; kdaMatrix=fkd; bNodeSz=40; bFontSz=12

      nodeprop = configureNodeVisualization(allnodes=allnodes, signature=signature, kdaMatrix=fkd)

      hnList     = nodeprop[[1]] # node subcategpries
      listprop   = nodeprop[[2]] # visual properties for each subcategory
      legend     = nodeprop[[3]] # legend table for visual propertie

      resf = makeSNP(netpairsWtype   = expandNet, 
               edgecolorlevels = c("grey"),
               highlightNodes  = hnList,
               normColor="grey",   highColor=listprop[,1],
               normShape="circle", highShape=listprop[,2],
               normNodeSize ="40",  highNodeSize =listprop[,3],
               normFontSize ="12",  highFontSize =listprop[,4],
               legendtable=legend, snafile=snpFname, kdaMatrix=nodeDegree)

      result = list(expandNet, fkd, ret[[2]], getFileFullNameNopath(resf) )
      names(result) <- c("subnetwork", "keydrivers", "parameters", "files")
   } else{
      result = list(expandNet, fkd, ret[[2]])
      names(result) <- c("subnetwork", "keydrivers", "parameters")
   }

   return (result)
}

keydriverInSubnetwork <- function(linkpairs, signature, background=NULL, directed=T, nlayers=6, 
                                   enrichedNodes_percent_cut=-1, FET_pvalue_cut=0.05, 
                                   boost_hubs=T, dynamic_search=T, bonferroni_correction=T, expanded_network_as_signature =F) {

   allnodes = sort(union(linkpairs[,1], linkpairs[,2]))
   no.subnetsize = length(allnodes)

   #print(no.subnetsize)
   if(no.subnetsize <=4) {return (NULL)}

   # whole network nodes as the signature
   network_as_signature = length(setdiff(allnodes, signature)) ==0
   network_as_signature = expanded_network_as_signature | network_as_signature

   overlapped    = intersect(allnodes, signature)
   no.overlapped = length(overlapped) # within the subnetwork

   if(is.null(background) ){
      background2 = c(no.subnetsize, no.overlapped) 
   } else {
      background2 = background
   }

   keydrivers= NULL
   kdMatrix  = NULL
   kdIndex   = NULL # indices of keydrivers in dsnodes_list

   dsnodes_list = as.list(rep(0,no.subnetsize)); no.dsnodes = rep(0, no.subnetsize)
   cnt = 1

   intv = as.integer(no.subnetsize/10)
   if(intv==0) {intv =1}
   print("find downstream genes")

   # set up searching range
   if(dynamic_search) { # dynamic search for optimal layer
      layers_4search = c(1:nlayers)
   } else{  # fixed layer for optimal layer
      layers_4search = c(nlayers)
   }
   # if the network itself is the signature, no need for dynamic search
   if(network_as_signature){layers_4search = c(nlayers)}

   dn_matrix = matrix(1, no.subnetsize, nlayers)
   for(i in c(1:no.subnetsize) ) {

     if(i%%intv==0){ print(paste(i, "/", no.subnetsize)) }

     # initialization
     minpv=1;min_nohits=0;min_noidn=0;min_layer=0;min_dn=0; min_fc=0;
     minpvW=1;min_fcW=0;
     for(y in layers_4search) {
         #netpairs=linkpairs; seednodes=allnodes[i]; N=nlayers; directed=directed
         idn = downStreamGenes(netpairs=linkpairs, seednodes=allnodes[i], N=y, directed=directed)
         idn = setdiff(idn, allnodes[i])
         no.idn=length(idn);

         dn_matrix[i, y] = no.idn
      
         if(!network_as_signature){# do enrichment test for only expanded subnetwork
            hits    = intersect(idn, overlapped)
            no.hits = length(hits)

            if(no.hits==0){next}

            foldchg = (no.hits/no.idn)/(no.overlapped/no.subnetsize)
            pv = phyper(no.hits-1, no.idn, no.subnetsize-no.idn, no.overlapped, lower.tail=F)

            foldchgW = (no.hits/no.idn)/(background2[2]/background2[1])
            pvW = phyper(no.hits-1, no.idn, background2[1]-no.idn, background2[2], lower.tail=F)

            if(pv<minpv){
               minpv=pv;min_nohits=no.hits;min_noidn=no.idn;min_layer=y;min_fc=foldchg; 
               minpvW=pvW;min_fcW=foldchgW
            }
         } else{ # for non-expanded subnetwork
            no.hits = no.idn
            minpv=0;min_nohits=no.idn;min_noidn=no.idn;min_layer=y;min_fc=1
         }  
     } #y
     
     # record the down stream genes for the biggest layer
     if(no.idn>0) {
        dsnodes_list[[i]] = idn
        no.dsnodes[i]     = no.idn
     }

     correct_minpv = minpv*no.subnetsize; correct_minpv = ifelse( correct_minpv>1, 1, correct_minpv)

     res= c(min_nohits,min_noidn,no.overlapped,no.subnetsize, background2[2], background2[1], 
            length(signature), min_layer,min_fcW, minpvW, min_fc,minpv, correct_minpv)
     kdMatrix = rbind(kdMatrix, res)
     #print(res)
   }

   colnames(kdMatrix) <- c("hits", "downstream", "signature_in_subnetwork", "subnetwork_size", 
                           "signature_in_network", "network_size", "signature", "optimal_layer",
                           "fold_change_whole", "pvalue_whole", "fold_change_subnet", "pvalue_subnet", 
                           "pvalue_corrected_subnet")
   optLidx = getMatchedIndexFast(colnames(kdMatrix), "optimal_layer")

   mymincut = enrichedNodes_percent_cut*no.overlapped
   if (enrichedNodes_percent_cut<=0){
      #mymincut = mean(no.dsnodes) + sd(no.dsnodes)
      mincutLayers = apply(dn_matrix, 2, mean) + apply(dn_matrix, 2, sd)
      opL = ifelse(kdMatrix[, optLidx]==0, 1, kdMatrix[, optLidx])
      mymincut = mincutLayers[opL]
   }
   cutmatrix = c( mean(no.dsnodes), sd(no.dsnodes), mymincut)

   # pick up key drivers by pvalue and no. of downstream genes
   ncols = dim(kdMatrix)[2]

   if(bonferroni_correction) { # use corrected pvalue
      kdSel = (kdMatrix[,ncols] < FET_pvalue_cut) & (kdMatrix[,2]>= mymincut)
   } else{
      kdSel = (kdMatrix[,ncols-1] < FET_pvalue_cut) & (kdMatrix[,2]>= mymincut)
   }
   
   keydrv = rep(0, no.subnetsize)
   if( sum(kdSel) >0){

      keydrivers= allnodes[kdSel]
      kdIndex   = c(1:no.subnetsize)[kdSel]
      n.drivers = length(keydrivers)

      #******************* local driver or not **************************************
      #
      # check whether a driver is in the downstream of other drivers
      #keydrv = rep(0, no.subnetsize)
      #if (!network_as_signature) {
      for ( i in c(1:n.drivers) ) {

          # Note that kdIndex[i] is the index of ith keydriver in kdMatrix  
          # restrict to only candidate drivers 
          iselA  = (kdMatrix[,2] > kdMatrix[ kdIndex[i],2]) & kdSel
          isel   = c(1:no.subnetsize)[iselA]

          if ( sum(isel)>0) {
              if(directed) {
                 ilocal= setInSets(setC=allnodes[ kdIndex[i] ],       setlist=dsnodes_list[isel])
              } else{
                 ilocal= setInSets(setC=dsnodes_list[[ kdIndex[i] ]], setlist=dsnodes_list[isel])
              }
              keydrv[ kdIndex[i] ] = !ilocal + 0
          } else{
              keydrv[ kdIndex[i] ] = TRUE
          }
      }
     #}
   } else{
      print("Warning: the downstream metric is not used as the specified minimumal downstream size is too big !!")
   }

   # promote genes with many direct links to be key drivers
   #
   #              inlinks outlinks totallinks
   #0610031J06Rik       2        0          2
   #1110001J03Rik       0        1          1
   #
   if(boost_hubs) {
  
     if(!network_as_signature){
        # for expanded network, restrict the boosted nodes to the key driver candidates
        kdSelB = rep(F, no.subnetsize); kdSelB[kdIndex]=T;
        
        psel = kdMatrix[,ncols-3]*no.subnetsize < 0.05; kdPsd=1; kdPmean =1;
        kdpvalues = -log10(kdMatrix[,ncols-3]); kdpvalues = ifelse(is.na(kdpvalues), 0, kdpvalues)
        #histogram(kdpvalues)

        if(sum(psel)>0) {
          kdPmean= mean(kdpvalues[psel]); kdPsd= sd(kdpvalues[psel])
          #kdPmean= median(kdpvalues[psel]); kdPsd= mad(kdpvalues[psel])

          print( as.numeric(signif(kdpvalues[psel],2)) )
          directSel = (kdpvalues > (kdPmean + kdPsd ) ); directSel=ifelse(is.na(directSel), FALSE, directSel)
          #print(directSel)
          if( sum(directSel)>0) {
            kdSel = kdSel | directSel
            dIndex    = c(1:no.subnetsize)[kdSel]
            keydrv    = rep(F, no.subnetsize)
            keydrv[dIndex] = TRUE
          }
        }
        cutmatrix = rbind( c(mean(no.dsnodes), sd(no.dsnodes), concatenate(mymincut,";"),kdPmean,kdPsd, kdPmean + kdPsd ))

        colnames(cutmatrix) <- c("mean_downstream", "sd_downstream", "enrichedNodes_cut", 
                                 "mean_logP", "sd_logP", "cut_logP")

     } else{
        # for non-expanded network, consider all the nodes in the subnetwork
        kdSelB = rep(TRUE, no.subnetsize);

        # align the degree with allnodes
        mydegree  = degreeByLinkPairs(linkpairs=linkpairs, directed=directed, cleangarbage=F)
        mIdx = getMatchedIndexFast(rownames(mydegree), allnodes)
        mydegree = mydegree[mIdx, ]

        if(directed) {
            directSel = mydegree[,2]> mean(mydegree[,2]) + 2*sd(mydegree[,2])
            cutmatrix = rbind( c(mean(no.dsnodes), sd(no.dsnodes), concatenate(mymincut, ";"),
                         mean(mydegree[,2]),sd(mydegree[,2]), mean(mydegree[,2]) + 2*sd(mydegree[,2]) ))
        }else{
            directSel = mydegree[,3]> mean(mydegree[,3]) + 2*sd(mydegree[,3])
            cutmatrix = rbind( c(mean(no.dsnodes), sd(no.dsnodes), concatenate(mymincut, ";"),
                        mean(mydegree[,3]),sd(mydegree[,3]), mean(mydegree[,3]) + 2*sd(mydegree[,3])))
        }
        directSel = directSel & kdSelB

        directeHub  = rownames(mydegree)[directSel]
        isDirectHub = setElementInSet(allnodes, directeHub)

        keydrv[isDirectHub] = T
        kdSel = kdSel | isDirectHub
        colnames(cutmatrix) <- c("mean_downstream", "sd_downstream", "cut_downstream",
                                 "mean_degree", "sd_degree", "cut_degree")
    }


   } else{
     cutmatrix = rbind( c(mean(no.dsnodes), sd(no.dsnodes), concatenate(mymincut, ";"), "F"))
     colnames(cutmatrix) <- c("mean_downstream", "sd_downstream", "cut_downstream", "boost_directhubs")
   }

   if( sum(kdSel)==0){return(NULL)}

   ##
   # in this case, signature is the network nodes themselves, so pvalue will be 0 for all nodes
   # so the driver will be the ones with most downsttream genes
   #
   is_signature = rep(0, no.subnetsize); names(is_signature) <- as.character(allnodes)
   is_signature[as.character(overlapped)]  = 1
   
   fkd = cbind(allnodes, is_signature, kdMatrix, keydrv+0)[kdSel,];

   if(sum(kdSel) >1 ) {
       nf.cols = dim(fkd)[2]
       if(network_as_signature){
          mo = order(-as.integer(fkd[,3]))
       } else{
          mo = order(as.numeric(fkd[,nf.cols-1]))
       }

       fkd = fkd[mo, ]
       # put key driver on the top
       mo  = order( -as.integer(fkd[,nf.cols]) )
       fkd = fkd[mo, ]
   } else{
       fkd = rbind(fkd)
   }

   colnames(fkd) <- c("keydrivers", "is_signature", "hits", "downstream", "signature_in_subnetwork", 
                      "subnetwork_size", "signature_in_network", "network_size", "signature", "optimal_layer",
                      "fold_change_whole", "pvalue_whole", "fold_change_subnet", "pvalue_subnet", "pvalue_corrected_subnet", "keydriver")

   print(fkd)

   degreeAll = cbind(allnodes, no.dsnodes); colnames(degreeAll)= c("node", "downstream")

   return( list( fkd , cutmatrix, degreeAll) )
}

makeSNP <- function( netpairsWtype , highlightNodes = NULL , edgecolorlevels , normColor = "grey" ,
		     highColor = "red" , normShape = "50" , highShape = "50", 
                     normNodeSize ="1",  highNodeSize="2", normFontSize ="1",    highFontSize="2",
                     directed = TRUE, legendtable = NA , snafile = "tmp.sna", kdaMatrix)
{
    xfname = getFileName(snafile)
    fcys  = paste(xfname, "_cys.txt", sep="")
    fcysn = paste(xfname, "_cys-nodes.txt", sep="")
    write.table(netpairsWtype, fcys, sep="\t", quote=FALSE, col.names=T, row.names=FALSE)

    pcols = dim(netpairsWtype)[2]

    # consider both columns
    uniquenames    = union(as.character(netpairsWtype[,1]), as.character(netpairsWtype[,2]) )
    uniquenames    = sort(uniquenames)
    no.uniquenames = length(uniquenames)

    # initialize matrix with the size equal to no.uniquenames
    #
    name2idxMatrix = cbind(uniquenames, c(1:no.uniquenames) )
    

    # 0. make link index: A1 A2 T & A I
    # A1 A2 T & A I ==> A1 A2 T I1
    #
    leftIdx = merge(netpairsWtype, name2idxMatrix, by.x=1, by.y=1, all=F) 
    
    #  A1 A2 T I1 & A I ==> A2 A1 T I1 I2: I1 and I2 are the indices of A1 and A2 respectively
    #
    allIdx  = merge(leftIdx, name2idxMatrix, by.x=2, by.y=1, all=F)

    no.pairs = dim(allIdx)[1]

    if(pcols==2){
      no.elevels = length(edgecolorlevels)
      greyIdx    = c(1:no.elevels) [edgecolorlevels=="grey"]
      linksIdxMatrix = cbind( allIdx[,c(3,4)], rep(greyIdx,no.pairs) )
    } else{
      linksIdxMatrix = allIdx[,c(4,5,3)]
    }


    # 1. head string
    header = paste("*vertices ", as.character(no.uniquenames), sep="")

    # 2. vertices matrix 
    if ( is.null(highlightNodes) ){
        verticesMatrix = cbind(uniquenames, 
                               rep(normColor,no.uniquenames), 
                               rep(normShape,no.uniquenames),
                               rep(normNodeSize, no.uniquenames),
                               rep(normFontSize, no.uniquenames))
    } else{
      verticesMatrix = matrix("", no.uniquenames, 5)
      verticesMatrix[,1] = uniquenames
      verticesMatrix[,2] = rep(normColor,no.uniquenames)
      verticesMatrix[,3] = rep(normShape,no.uniquenames)
      verticesMatrix[,4] = rep(normNodeSize, no.uniquenames)
      verticesMatrix[,5] = rep(normFontSize, no.uniquenames)

      xcolor= rep(normColor,no.uniquenames);    names(xcolor) <- uniquenames
      xshape= rep(normShape,no.uniquenames);    names(xshape) <- uniquenames
      xNsize= rep(normNodeSize,no.uniquenames); names(xNsize) <- uniquenames
      xSsize= rep(normFontSize,no.uniquenames); names(xSsize) <- uniquenames

      # set color and shape for highlighted nodes
      #
      if ( !is.list(highlightNodes) ){
          highlightNodes2 = intersect(highlightNodes, uniquenames)
          xcolor[ highlightNodes2] = highColor[1]
          xshape[ highlightNodes2] = highShape[1]
          xNsize[ highlightNodes2] = highNodeSize[1]
          xSsize[ highlightNodes2] = highFontSize[1]
      }else{
         no.highnodeSets = length(highlightNodes)
         for(il in c(1:no.highnodeSets) ) {
             highlightNodes2 = intersect(highlightNodes[[il]], uniquenames)
             if(length(highlightNodes2)==0){next};
             xcolor[highlightNodes2] = highColor[il]
             xshape[highlightNodes2] = highShape[il]
             xNsize[highlightNodes2] = highNodeSize[il]
             xSsize[highlightNodes2] = highFontSize[il]
         }
      }
      verticesMatrix[,2] = as.character(xcolor)
      verticesMatrix[,3] = as.character(xshape)
      verticesMatrix[,4] = as.character(xNsize)
      verticesMatrix[,5] = as.character(xSsize)
    }

    #verticesMatrix <- as.matrix(verticesMatrix)
    colnames(verticesMatrix) <- c("nodename","color","shape", "size", "font_size")

    verticesMatrix2= verticesMatrix[,c(1:3)]
    verticesMatrix2= merge(verticesMatrix2, kdaMatrix, by.x=1, by.y=1, all.x=TRUE, sort=FALSE)
    verticesMatrix3= as.matrix(verticesMatrix2)[,c(1:ncol(verticesMatrix2), ncol(verticesMatrix2)) ]
    colnames(verticesMatrix3) <- c("nodename","color","shape", "size", "font_size")
    write.table(verticesMatrix3, fcysn, sep="\t",quote=FALSE, col.names=T, row.names=FALSE)

    #**************************************************************************
    #
    # 3. output indexed netpairs
    #
    # Legend
    if ( !is.na(legendtable) ) {
      mhead = paste("Legend", dim(legendtable)[1], sep=" ")
      write.table(as.matrix(mhead), snafile, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=F)
      write.table(legendtable, snafile, sep="\t",quote=FALSE, col.names=T, row.names=FALSE, append=T)
      # vertex
      write.table(as.matrix(header), snafile, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T)
    } else{
      # vertex
      write.table(as.matrix(header), snafile, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=F)
    }

    write.table(verticesMatrix, snafile, sep="\t",quote=FALSE, col.names=T, row.names=FALSE, append=T)

    #edge color
    write.table(t(as.matrix(c("edge_color_levels", edgecolorlevels) )), snafile, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T)

    # link pairs based index   
    #
    write.table(rbind(c("src","dst", "type")), snafile, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T)
    write.table(linksIdxMatrix, snafile, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T)

    return (c(fcys, fcysn))
}

mergeTwoMatricesByKeepAllPrimary <- function( primaryMatrix , minorMatrix , missinglabel = "" ,
		                            keepAllPrimary = TRUE , keepPrimaryOrder = TRUE ,
									keepAll = FALSE )
{
  no.promarycols <- dim( primaryMatrix )[2]
  no.mustbegenes <- dim( primaryMatrix )[1]

  # we add in one more column to indicate which genes are mustbeincluded after being merged with mcg
  keyword <- "mustbeused"
  mustbeGenesMatrix <- cbind( primaryMatrix , c( 1:no.mustbegenes ) , rep( keyword , no.mustbegenes ) )

  if ( is.null( colnames( primaryMatrix ) ) )
  {
    colnames( mustbeGenesMatrix ) <- c( c( 1:no.promarycols ) , "primorder" , keyword )
  }
  else
  {
    colnames( mustbeGenesMatrix ) <- c( colnames( primaryMatrix ) , "primorder" , keyword )
  }
# Why is this uncommented?
#  dim( mustbeGenesMatrix )

  if ( is.null( keepAllPrimary ) )
  { #normal merge: to have the common elements
    myMatrix <- merge( mustbeGenesMatrix , minorMatrix , by.x = 1 , by.y = 1 , all.x = FALSE ,
			           sort = FALSE , all = FALSE )
  }
  else
  {
    myMatrix <- merge( mustbeGenesMatrix , minorMatrix , by.x = 1 , by.y = 1 , all.x = TRUE ,
			           sort = FALSE , all = TRUE )
  }
# Again, why is this left uncommented?
#  dim( myMatrix )
  nocols.mymatrix <- dim( myMatrix )[2]

  #the mustbeused genes which are not included in minor have NAs in the column $mustbeused
  #so we can use this information to figure out which mustbeused genes missing in minorMatrix
  myMatrix[,nocols.mymatrix] <- ifelse( is.na( myMatrix[,nocols.mymatrix] ) , missinglabel ,
		                        as.character( myMatrix[,nocols.mymatrix] ) )

  orders <- order( as.numeric( as.matrix( myMatrix[,no.promarycols + 1] ) ) )
  if ( keepPrimaryOrder )
      myMatrix <- myMatrix[orders,]

  if ( is.null( keepAllPrimary ) )
  {
     selected <- rep( T , dim( myMatrix )[1] )
  }
  else
  {
     if ( keepAllPrimary )
	 {
       selected <- !( is.na( myMatrix[,no.promarycols + 2] ) )
     }
     else #return the row-elements in minor which are missed in primary
	 {
	   selected <- is.na( myMatrix[,no.promarycols + 2] )
	 }
  }
# Why are these here?
#  sum(selected)

  #keep the primary matrix and remove the mustbeused column
  return( myMatrix[selected,-c( no.promarycols + 1 , no.promarycols + 2 )] )
}

removeDuplicatedLinks <- function( linkpairs , directed = FALSE )
{
    if ( isTRUE( all.equal( dim( linkpairs )[1] , 1 ) ) )
	{
       return( linkpairs )
    }

    links <- paste( linkpairs[,1] , linkpairs[,2] , sep = "\t" )

    # 1. remove duplications 
    #
    cleanedlinkMatrix <- union( links , NULL )
# Why do we keep having these random printing things??
#	length( cleanedlinkMatrix )

    #ofname ="tmp.txt"
    #write.table( cleanedlinkMatrix, ofname, sep="\t", quote=FALSE, col.names=F, row.names=FALSE)

    # 2. remove inversed duplications
    #
    #linkMatrix <- read.delim(ofname, sep="\t", header=F) # not good for batch operation, ie, many of the same jobs running
    linkMatrix  <- getAllParts( cleanedlinkMatrix , "\t" )
# Again, another one!
#    dim(linkMatrix)
    linkMatrix <- as.matrix( linkMatrix )

    if ( directed )
	{
       return( linkMatrix )
    }

    if ( isTRUE( all.equal( dim( linkMatrix )[1] , 1 ) ) )
	{
		return( linkMatrix )
	}

    #  first, remove self-interactions are also removed
    #
    selSelfLinks <- linkMatrix[,1] == linkMatrix[,2]
    linkMatrix <- linkMatrix[!selSelfLinks,]
    cleanedlinkMatrix <- cleanedlinkMatrix[!selSelfLinks]

    reversedLinks <- paste( linkMatrix[,2] , linkMatrix[,1] , sep = "\t" )

    no.links <- length( reversedLinks )
    reversedLinksWithIndex <- cbind( reversedLinks , c( 1:no.links ) )

    merged <- merge( cleanedlinkMatrix , reversedLinksWithIndex , by.x = 1 , by.y = 1 , all = FALSE )
# What is it with these things?!?
#	dim(merged)

    if ( dim( merged )[1] > 0 )
	{
        merged <- as.matrix( merged )
        removedCols <- as.integer( merged[,2] )

        # cosntruct non-duplicated interactions
        #
        dupLinks <- cleanedlinkMatrix[removedCols]
        dupLinksRev <- reversedLinksWithIndex[removedCols]
        uniques <- NULL
        for ( i in c( 1:length( dupLinks ) ) )
        {
           found <- is.element( dupLinks[i] , uniques )
           foundRev <- is.element( dupLinksRev[i] , uniques )
           combined <- found | foundRev
           if ( !combined )
		   {
               uniques <- c( uniques , dupLinks[i] )
           }
        }
# ANOTHER ONE?!?
#		length(uniques)
        xlinkMatrix <- c( cleanedlinkMatrix[-removedCols] , uniques )
        #write.table( cleanedlinkMatrix[-removedCols], ofname, sep="\t", quote=FALSE, col.names=F, row.names=FALSE)
        #write.table( uniques, ofname, sep="\t", quote=FALSE, col.names=F, row.names=FALSE,append=T)
    }
	else
	{
        #write.table( cleanedlinkMatrix, ofname, sep="\t", quote=FALSE, col.names=F, row.names=FALSE)
        xlinkMatrix <- cleanedlinkMatrix
    }

    #linkMatrix <- read.delim(ofname, sep="\t", header=F)
    #dim(linkMatrix)
    #linkMatrix  = as.matrix(linkMatrix)

    return( getAllParts( xlinkMatrix , "\t" ) )
}

replaceString <- function( fullfnames , oldstr , newstr )
{
  no.files <- length( fullfnames )
  res <- NULL
  for ( each in fullfnames )
  {
    #print(paste(i, "/", no.files, ":", each) )
    each2 <- paste( each , oldstr , sep = "" )
    splitted <- splitString( each2 , oldstr )

    neweach <- concatenate( splitted , newstr )

# How could this ever be run?  There is no way for this conditional to evaluate to TRUE!
#    if ( FALSE )
#	{
#      neweach <- ""
#      for ( is in splitted )
#	  {
#        neweach <- paste( neweach , newstr , sep = is )
#      }
#    }

    #oldeach  = paste(pathnet, each,   sep="")
    #neweach  = paste(pathnet, newstr, splitted[2], sep="")
    #a=file.rename(from=oldeach, to=neweach)
    #print(a)

    res <- c( res , neweach )
  }
  return( res )
}

setElementInSet <- function( setA , setB )
{
    found <- rep( FALSE , length( setA ) )
    for ( i in c( 1:length( setA ) ) )
	{
       idiff <- setdiff( setA[i] , setB )
	   found[i] <- isTRUE( all.equal( length( idiff ) , 0 ) )
    }
    return ( found )
}

setInSets <- function( setC , setlist )
{
   for ( i in c( 1:length( setlist ) ) )
   {
       isSub <- setsub( setC , setlist[[i]] )
       if ( isSub )
	   { #setC is a subset of setlist[[i]]
          return( TRUE )
       }
   }
   return( FALSE )
}

setsub <- function( setA , setB )
{
    if ( length( setA ) > length( setB ) )
	{
        return( FALSE )
    }
    setAB <- union( setA , setB )
    return ( setequal( setAB , setB ) )
}

splitString <- function( mystring , separator = "; " )
{
  splitted <- NULL
  for ( each in mystring )
  {
     if ( is.na( each ) | is.null( each ) )
	 {
        next
     }
     splitted <- c(splitted , unlist( strsplit( each , separator ) ) )
  }
  #a=unlist( strsplit(mystring, separator) )
  return( splitted )
}

