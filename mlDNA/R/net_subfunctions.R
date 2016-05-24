###########################################################################
#####threshold is used to remove low correlation from adjacy matrix

######################################################
##matrix permutation
.permutMatrix <- function( mat, permutByRow = NULL, permutByCol = NULL ) {
  if( is.null(permutByRow) & is.null(permutByCol) ){
    stop("Error: no permutation?")
  }
  
  .permutVec <- function(vec, len){
    idx <- sample(len, replace=FALSE)
    return(vec[idx])
  }
  
  nrows <- nrow(mat)
  ncols <- ncol(mat)
  
  newmat <- NULL
  
  if( !is.null(permutByRow) & is.null(permutByCol) )
    newmat <- t(apply( mat, 1, .permutVec, len=ncols))
  if( !is.null(permutByCol) & is.null(permutByRow) )
    newmat <- apply( mat, 2, .permutVec, len=nrows)
  
  if( !is.null(permutByRow) & !is.null(permutByCol) ){
    newmat <- t(apply( mat, 1, .permutVec, len=ncols))
    newmat <- apply( newmat, 2, .permutVec, len=nrows)
  }
  
  newmat
  
}


.permutThreshold <- function( permutVec, sigLevels = c(0.05, 0.01), side.type = c("one.sided", "two.sided" ) ) {
  
  vec <- NULL
  if( side.type == "one.sided") {
    vec <- permutVec
  }else {
    vec <- abs(permutVec)
  }
  vec <- sort(permutVec, decreasing=TRUE)
  
  len <- length(vec)
  idx <- ceiling(len*sigLevels)
  if( length(sigLevels) == 1 ){
    if( idx <= 0 )
      idx <- 1
  }else{
    idx[idx < 0] <- 1
  }
  t <- vec[idx]
  names(t) <- as.character(sigLevels)
  t
}




##################################################################################
##SubFunction: get threshold of correlation coefficient with the method presented in 
##Bioinformatics, 2004, 20(14): 2242-2250.
##################################################################################
.cor.threshold <- function( expMat, sigLevels = c(0.05, 0.01, 0.001), corMethod = c("GCC", "PCC", "SCC", "KCC", "BiWt", "MI", "MINE"), 
                            distMethod = c("tDist", "permutation"), tailed = "two.sided", cpus = cpus, 
                            saveType = "matrix", backingpath= NULL, backingfile= NULL, descriptorfile= NULL  ) {
  
  samplesize = ncol(expMat)
  res <- rep(0, length(sigLevels))
  names(res) <- sigLevels
  
  method <- distMethod
  if( is.null(distMethod) )
    stop("Error: distMethod should be defined")
  if( is.vector(distMethod) )
    method <- distMethod[1]
  
  if( method == "tDist" ){
    for( i in 1:length(sigLevels) ){
      tTh <- qt(1-sigLevels[i], samplesize - 2, lower.tail = TRUE, log.p = FALSE)
      res[i] <- 1.0*tTh/sqrt(samplesize-2+tTh^2)
    }
  }else{
    if( is.null(corMethod) )
      stop("Error: corMethod should be specified when distMethod is \"permutation\" ")
    
    expMat.permut <- .permutMatrix( expMat, permutByRow = TRUE, permutByCol = NULL )  #y independently permuting the components of each gene expression vector
    permutMat <- adjacencymatrix( expMat.permut, method= corMethod, cpus = cpus, saveType = saveType, backingpath= backingpath, backingfile= backingfile, descriptorfile= descriptorfile )
    permutMat1 <- permutMat[,]
    permutVec <- sort( permutMat1[upper.tri(permutMat1)], decreasing= TRUE ) 
    
    res <- .permutThreshold( permutVec, sigLevels =sigLevels, side.type = tailed )
    remove("permutMat")
    remove("expMat.permut")
    remove("permutVec")
  }
  res
}


#################################################
##filter low correlation
.filterLowCorrelation <- function( adjmat, threshold ) {
  t <- NULL
  if( is.big.matrix(adjmat) ) {
    t <- adjmat[,]
  }else{
    t <- adjmat
  }
  t[ which(t > (-1.0*threshold) & t < threshold)] <- 0
  t
}



##########################################################################################
##SubFun: format adjmatrix to graph
#######################################################################################
.formatAdj2Graph <- function( adjmat, threshold, file = NULL, format = c("edgelist", "pajek", "ncol", "lgl", "graphml", "dimacs", "gml", "dot", "leda"), ... ) {
  
  call <- match.call()
  
  if( length(format) > 1 & !is.null(file)) 
    cat("Warnning: only the first one format will be saved.")
  format <- format[1]
  
  t <- .filterLowCorrelation( adjmat = adjmat, threshold = threshold)
  igraph <- graph.adjacency( abs(t), mode="max", weighted= TRUE, diag=FALSE)
  
  
  if( !is.null(file) ) {
    write.graph( graph = igraph, file = file, format = format, ...)
  }
  
  rm(t)
  
  igraph
}


#get indices of connected genes, 20130616
.DNA.connectivitylist <- function( adjmat, threshold, backingpath = NULL, descriptorfile = NULL, nodes =NULL, cpus = 1 ) {
  
  if( is.null( rownames(adjmat)) )
    stop("Error: no rownames of adjmat")
  
  descfile <- .checkadjmatrix(mat = adjmat, backingpath = backingpath, descriptorfile = descriptorfile)
  nodes.idx <- NULL
  if( !is.null(nodes) ) {
    nodes.idx <- .getIndex( rownames(adjmat), nodes )
  }else{
    nodes.idx <- 1:nrow(adjmat)
  }
  
  
  .fun_conn_idx <- function( nodeIndex, mat, threshold, descfile ) {
    v <- NULL
    if( !is.null(descfile) ){
      v <- attach.big.matrix(descfile)[nodeIndex,]
    }else {
      v <- mat[nodeIndex,]
    }
    #mask diag
    v[nodeIndex] <- 0 
    pos.idx <- which(v >= threshold)
    neg.idx <- which(v <= -1.0*threshold)
    all.idx <- c(pos.idx, neg.idx)
    return( list(pos = pos.idx, neg = neg.idx, all = all.idx ) )
  }
  
  if( cpus == 1 ) {
    result <- apply( matrix(nodes.idx, ncol=1), 1, .fun_conn_idx, mat = adjmat, threshold = threshold, descfile = descfile)
  }else {
    sfInit(parallel = TRUE, cpus = cpus)
    sfLibrary( "bigmemory", character.only=TRUE )
    result <- sfApply( matrix(nodes.idx, ncol=1), 1, .fun_conn_idx, mat = adjmat, threshold = threshold, descfile = descfile)
    sfStop() 
  }
  
  names(result) <- nodes
  result
  
}



###################################################################################################
########network properties########################
.network.properties <- function(nodes = NULL, adjmat, threshold,  backingpath = NULL, backingfile = NULL, descriptorfile = NULL,  
                               graph = NULL, distmat = NULL, connectivityList = NULL, knodes = NULL, cpus = 1, 
                               properties = c("AllConnectivity", "PosConnectivity", "NegConnectivity", "1N", "2N",
                                              "closeness", "eccentricity", "eigenvector", "page.rank", 
                                              "dis2knodes", "closeness2knodes", "eccenticity2knodes"),
                               netDescribe = NULL, verbose = TRUE ) {
  
 
  properties_org <- properties
  propmat <- matrix(0, nrow = length(nodes), ncol = 1)
  rownames(propmat) <- nodes
  colnames(propmat) <- "tmp"
  ###########################################################################
  ##for check connectivities
  connProp <- c("AllConnectivity", "PosConnectivity", "NegConnectivity")
  overlapProp <- intersect( connProp, properties )
  if( length(overlapProp) >= 1 ) {
    if( verbose )
      cat("\n...start to calcluate properties: ", overlapProp, "for netowrk ", netDescribe, "...\n")
    if( is.null(connectivityList) ) {
        connectivityList <- .DNA.connectivitylist( adjmat = adjmat, threshold = threshold, backingpath = backingpath, descriptorfile = descriptorfile, nodes = nodes, cpus = cpus )  
        save( connectivityList, file = paste( backingpath, netDescribe, "_graph_connectivityList.RData", sep = "") )
    }
    conRes <- .DNA.ConnectivityNum( nodes = nodes, ConnectivityList =  connectivityList, cpus = cpus )
    conRes <- conRes[,overlapProp]
    if( !is.matrix(conRes) ) {
      conRes <- matrix( conRes, ncol = 1 )
      rownames(conRes) <- nodes
      colnames(conRes) <- overlapProp
    }
    curColNames <- c( colnames(propmat), overlapProp )
    propmat <- cbind( propmat, conRes[,overlapProp] )
    colnames(propmat) <- curColNames
    properties <- setdiff( properties, overlapProp )
  }

  ###########################################################################
  connProp <- c("1N", "2N")
  overlapProp <- intersect( connProp, properties )
  if( length(overlapProp) >= 1 ) {
    if( verbose )
      cat("\n...start to calcluate properties: ", overlapProp, "for network ", netDescribe, "...\n")
    if( is.null(connectivityList) ) {
      connectivityList <- .DNA.connectivitylist( adjmat = adjmat, threshold = threshold, backingpath = backingpath, descriptorfile = descriptorfile, nodes = nodes, cpus = cpus )  
      save( connectivityList, file = paste( backingpath, netDescribe, "_graph_connectivityList.RData", sep = "") )
    }
    tmpMat <- .DNA.NIndex2knodes( adjmat = adjmat, connectivitylist = connectivityList, nodes = nodes, knodes = knodes, cpus = cpus, types = overlapProp )
    curColNames <- c( colnames(propmat), overlapProp )
    propmat <- cbind( propmat, tmpMat ) 
    colnames(propmat) <- curColNames 
    properties <- setdiff( properties, overlapProp )
  }
  
     
  ###########################################################################
  ##for eigenvector and page.rank
  connProp <- c( "eigenvector", "page.rank" )
  overlapProp <- intersect( connProp, properties )
  if( length(overlapProp) >= 1 ) {
    if( verbose )
      cat("\n...start to calcluate properties: ", overlapProp, "for network ", netDescribe, "...\n")
    tmp <- .DNA.getNodeProperty( genes = nodes, graph = graph, types= overlapProp, normalized = FALSE ) 
 
     curColNames <- c( colnames(propmat), overlapProp )
     propmat <- cbind( propmat, tmp ) 
     colnames(propmat) <- curColNames 
     properties <- setdiff( properties, overlapProp )
  }
  
  
  
  
  ###########################################################################
  #new properties.
  connProp <- c("closeness", "eccentricity", "dis2knodes", "closeness2knodes", "eccenticity2knodes")
  overlapProp <- intersect( connProp, properties )
  if( length(overlapProp) >= 1 ) {
    if( verbose )
      cat("\n...start to calcluate properties: ", overlapProp, "...\n")
    tmpMat <- matrix(NA, nrow = length(nodes), ncol = length(overlapProp))
    rownames(tmpMat) <- nodes
    colnames(tmpMat) <- overlapProp
        
    ##check whether some nodes have no connection
    AllNodesInGraph <- V(graph)$name
    if( nrow(adjmat) != length(AllNodesInGraph) )
      cat( nrow(adjmat) - length(AllNodesInGraph), " nodes have no connection, and thus are not included in network.")
    
    newnodes <- intersect( nodes, AllNodesInGraph)
    newknodes <-intersect( knodes, AllNodesInGraph)
    #generate distance matrix
    if( is.null(distmat) ) {
       backingfile_dist <- paste( netDescribe, "_graph_distmat_bfile", sep = "")
       descriptorfile_dist <- paste( netDescribe, "_graph_distmat_dfile", sep = "")
       system.time( distmat <- .DNA.distance( graph = graph, v = newnodes, to = AllNodesInGraph, cpus = cpus, saveType = "bigmatrix", backingpath = backingpath, backingfile = backingfile_dist, descriptorfile = descriptorfile_dist ) )
    }
    
    .checkDistmat( distmat = distmat, nodes = newnodes, knodes = newknodes )
   
    
    #for different parameters
    for( curProp in overlapProp ) {
      if( verbose )
         cat( "Current...", curProp, "\n")
      if( curProp == "closeness" ) {
        tmpMat[newnodes, curProp] <- .DNA.closeness( subDistMat = distmat[newnodes,], cpus = cpus )
      }else if( curProp == "eccentricity"  ) {
        tmpMat[newnodes,curProp] <- .DNA.eccentricity( subDistMat = distmat[newnodes,], cpus = cpus )
      }else if( curProp == "dis2knodes" ){
        tmpMat[newnodes,curProp] <- .DNA.distance2knowngenes( subDistMat = distmat[newnodes, newknodes], mean = TRUE, cpus = cpus )
      }else if( curProp == "closeness2knodes" ) {
        tmpMat[newnodes, curProp] <- .DNA.closeness( subDistMat = distmat[newnodes, newknodes], cpus = cpus )
      }else if(  curProp == "eccenticity2knodes" ) {
        tmpMat[newnodes,curProp] <- .DNA.eccentricity( subDistMat = distmat[newnodes, newknodes], cpus = cpus )
      }else {
        stop("Error, undefined network property")
      }   
    }#end for
  
    propmat <- cbind( propmat, tmpMat ) 
    properties <- setdiff( properties, overlapProp )
  }
  
  ###########################################################################
  if( length(properties) > 1 ) {
    cat(properties, " --> these propertes are not defined.\n")
  }
  
  propmat <- propmat[,-1]
  finalProperties <- setdiff( properties_org, properties )
  propmat <- propmat[, finalProperties]
  rownames(propmat) <- nodes
  propmat[which(propmat == NaN)] <- NA
   
  propmat
  
}




#############################################################################################
##features from the differences between two networks
.network.differences <- function( expmat1, net1, threshold1, backingpath1 = NULL, descriptorfile1 = NULL, connectivityList1 = NULL,  
                                 expmat2, net2, threshold2, backingpath2 = NULL, descriptorfile2 = NULL, connectivityList2 = NULL,
                                 nodes, properties = c("expDistance", "corDistance", "ASC"), cpus = 1, verbose = TRUE ) {
  
  descfile1 <- .checkadjmatrix(net1, backingpath1, descriptorfile1)
  descfile2 <- .checkadjmatrix(net2, backingpath2, descriptorfile2)
  
  if( is.null(connectivityList1) ){
    connectivityList1 <- .DNA.connectivitylist( adjmat = net1, threshold = threshold1, backingpath = backingpath1, descriptorfile = descriptorfile1, nodes = nodes, cpus = cpus )  
    save( connectivityList1, file = paste( backingpath1, descriptorfile1, "_connectivity.RData", sep = "") )
  }
  if( is.null(connectivityList2) ){
    connectivityList2 <- .DNA.connectivitylist( adjmat = net2, threshold = threshold2, backingpath = backingpath2, descriptorfile = descriptorfile2, nodes = nodes, cpus = cpus )  
    save( connectivityList2, file = paste( backingpath2, descriptorfile2, "_connectivity.RData", sep = "") )
  }
  
  properties_org <- properties
   propmat <- matrix(0, nrow = length(nodes), ncol = 1)
  rownames(propmat) <- nodes
  colnames(propmat) <- "tmp"
  
  ####################################################################
  ##for difExp
  if( is.element("expDistance", properties) ) {
    if(verbose) {
      cat("\n...start to calcluate properties: expDistance...\n") 
    }
    tmp <- (expmat1[nodes,] - expmat2[nodes,])^2
    tmp <- matrix( sqrt( apply( tmp, 1, mean ) ), ncol = 1 )
 
    curColNames <- c( colnames(propmat), "expDistance" )
    propmat <- cbind( propmat, tmp)
    colnames(propmat) <- curColNames
    properties <- setdiff( properties, "expDistance")
  }
  if( length(properties) == 0 ) {
    return(propmat)
  }

   if(verbose) {
     cat("\n...start to calcluate properties: ", properties, "...\n") 
   }
  

  #sub function
  .fun_distance <- function( node, expmat1, net1, descfile1, connectivityList1, 
                                   expmat2, net2, descfile2, connectivityList2) {

    res <- rep(0, 2)
    names(res) <- c("ASC", "corDistance")
    
    v1 <- names(connectivityList1[[node]]$all)
    v2 <- names(connectivityList2[[node]]$all)
    vv <- union( v1, v2 )
    res["ASC"] <- (2*length(vv) - length(v1) - length(v2))/2
    
    ###for other method    
    corVec1 <- NULL
    corVec2 <- NULL
    if( !is.null(descfile1) & !is.null(descfile2) ){
      corVec1 <- attach.big.matrix(descfile1)[node, vv]
      corVec2 <- attach.big.matrix(descfile2)[node, vv]
    }else{
      corVec1 <- net1[node, vv]
      corVec2 <- net2[node, vv]
    }
    dw_power2 <- (corVec1 - corVec2)^2
    
    ##for corDistance
    res["corDistance"] <- sqrt( mean(dw_power2) )
    res 
  }#end for fun
  
  
 
  result <- NULL
  if( cpus == 1 ) {
    system.time( result <- apply( matrix(nodes, ncol=1), 1, .fun_distance, expmat1 = expmat1, net1 = net1, descfile1 = descfile1, connectivityList1 = connectivityList1, 
                                                                           expmat2 = expmat2,  net2 = net2, descfile2 = descfile2, connectivityList2 = connectivityList2 ) )
  }else {
    sfInit(parallel = TRUE, cpus = cpus)
    sfLibrary("bigmemory", character.only=TRUE)
    system.time( result <- sfApply( matrix(nodes, ncol=1), 1, .fun_distance, expmat1 = expmat1,  net1 = net1, descfile1 = descfile1, connectivityList1 = connectivityList1, expmat2 = expmat2,  
                                                                             net2 = net2, descfile2 = descfile2, connectivityList2 = connectivityList2 ) )
    sfStop()
  }
  
  result <- t(result)
  rownames(result) <- nodes
 

  curColNames <- c( colnames(propmat), properties )
  propmat <- cbind( propmat, result[, properties])
  colnames(propmat) <- curColNames

  propmat[,-1]
}





##################################################################################
##SubFunction: get structural propterties of centralities
##################################################################################
.getNodeProperty <- function( genes, graph, 
                              types=c("alpha.centrality","betweenness", "closeness", "degree", "eccentricity", "eigenvector","page.rank", "subgraph.centrality"), 
                              normalized = TRUE) {
  
  m <- matrix(0, nrow = length(genes), ncol = length(types) )
  colnames(m) <- types
  rownames(m) <- genes
  
  
  for( ii in seq(length(types)) ){
    
    tt <- NULL
    if( types[ii] == "alpha.centrality") {
      cat("Warnning: function may be interrupt for the sigular computation for alpha.centrality parameter.\n")
      tt <- alpha.centrality(graph, genes, alpha=1, loops=FALSE, exo=1, weights=NULL, tol=1e-7, sparse=TRUE)
    }else if( types[ii] == "betweenness" ) {
      tt <- betweenness(graph=graph, v=genes, directed = FALSE, nobigint = TRUE, normalized = normalized)
    }else if( types[ii] == "closeness" ) { 
      tt <- closeness(graph=graph, vids = genes, mode = "all", normalized = normalized)
    }else if( types[ii] == "degree" ) {
      tt <- degree(graph=graph, v = genes, mode = "all", normalized = normalized)
    }else if( types[ii] == "eccentricity" ) {
      tt <- eccentricity(graph=graph, vids = genes, mode= "all")
    }else if( types[ii] == "page.rank" ) {
      tt <- page.rank (graph, vids = genes, damping = 0.85)$vector
    }else if( types[ii] == "eigenvector"  ) {
      tt <- evcent(graph, scale = normalized)$vector[genes]
    }else if( types[ii] == "subgraph.centrality" ){
      cat("Warnning: Much time would be took for huge network.\n")
      tt <- subgraph.centrality(graph, diag=FALSE)[genes]
    }else{
      stop("Undefined types in function:getCentrality")
    }
    
    if( length( which( (as.character(names(tt)) == as.character(genes)) == FALSE ) ) > 0 ) {
      cat(length(tt), "VS", length(genes), '\n')
      stop("Conflict between the order of gene names.\n")
    }
    
    m[1:length(genes),ii] <- tt
    
    
  }#end for 
  return (m)  
}






##get connectivity number
.DNA.ConnectivityNum <- function( nodes = NULL,  ConnectivityList, cpus = 1 ) {
  
  if( is.null(nodes) ) 
    nodes <- names(ConnectivityList)
  
  .subConnection <- function( node, ConnectivityList ) {
    tt <- ConnectivityList[[node]]
    vec <- rep(0, 3)
    vec[1] <- length( tt$all )
    vec[2] <- length( tt$pos )
    vec[3] <- length( tt$neg )
    vec
  }
  
  res <- NULL
  if( cpus  == 1 ) {
    res <- apply( matrix( nodes, ncol = 1), 1, .subConnection, ConnectivityList = ConnectivityList)
  }else {
    sfInit(parallel = TRUE, cpus = cpus)
    res <- sfApply( matrix( nodes, ncol = 1), 1, .subConnection, ConnectivityList = ConnectivityList )
    sfStop()   
  }
  res <- t(res)
  rownames(res) <- nodes
  colnames(res) <- c("AllConnectivity", "PosConnectivity", "NegConnectivity")
  res
}







##################################################################################
##SubFunction: get structural propterties of centralities
#note: use cpus= 1 for eigenvector and page.rank
##################################################################################
.DNA.getNodeProperty <- function( graph, genes, 
                             types=c("alpha.centrality", "betweenness", "closeness", "degree", "eccentricity", "page.rank", "eigenvector", "subgraph.centrality"), 
                             normalized = TRUE, cpus = 1) {
  
  if (cpus == 1 ) {
    result <- .getNodeProperty( graph = graph, genes = genes, types = types, normalized = normalized )
  }
  else {
    #if( !require(snowfall) ){
    #  install.packages("snowfall")
    #  library(snowfall)
    #}
    
    #not compute "eigenvector", "subgraph.centrality"
    types.new <- types[!( (types == "eigenvector") | (types == "subgraph.centrality") )]
    
    if( length(types.new) > 0 ) {
      genelist.matrix <- matrix(genes, ncol=1 )
      sfInit(parallel = TRUE, cpus = cpus)
      sfLibrary("igraph", character.only=TRUE )
      results.new <- sfApply( genelist.matrix, 1, .getNodeProperty, graph = graph,  types = types.new, normalized = normalized )
      sfStop()
    }
    
    result <- matrix(0, nrow=length(types), ncol=length(genes) )
    tt <- NULL
    for( ii in c(1:length(types))) {
      curType = types[ii]
      if( (curType == "eigenvector") | (curType == "subgraph.centrality")  ) {
        tt <- .getNodeProperty( graph = graph, genes = genes, types = curType, normalized = normalized )
      }else{
        tt <- results.new[which(types.new == curType),]
      }
      result[ii,] <- tt
    }
    result <- t(result)
  }
  
  if( length( types ) == 1 ){
    results <- matrix(result, ncol=1)
  }
  rownames(result) <- genes
  colnames(result) <- types
  result
}



##get distance between two genes
.DNA.distance <- function( graph, v, to, cpus = 1, saveType = "bigmatrix", backingpath = NULL, backingfile = NULL, descriptorfile = NULL){
  
  #if( !require(igraph) )
  #  library(igraph)
  
  if( saveType == "bigmatrix" ) {
    if( is.null(backingpath) | is.null(backingfile) | is.null(descriptorfile) )
      stop("Error: backingpath, backingfile and descriptorfile MUST be specified for bigmatrix")
  }
  
  .sub_DNA.distance <- function( node, graph, to ) {
    res <- shortest.paths( graph = graph, v = node, to = to, mode = "all", weights = NULL )
    res
  }
  
  ##calculate distance with weighted edges
  if( cpus == 1 ) {
    distmat = shortest.paths( graph, v = v, to = to, mode = "all", weights = NULL)  
  }else {
    sfInit(parallel = TRUE, cpus = cpus)
    sfLibrary("igraph", character.only=TRUE)
    distmat <- t( sfApply( matrix(v, ncol = 1), 1, .sub_DNA.distance, graph = graph, to = to  ) )
    sfStop()
  }
  diag(distmat) <- NA
  distmat[which(distmat == Inf)] <- NA
  distmat[which(distmat == -Inf)] <- NA
  rownames(distmat) <- v
  colnames(distmat) <- to
  
  if( saveType == "bigmatrix" )
    distmat <- as.big.matrix( distmat, backingfile = backingfile, backingpath = backingpath, descriptorfile = descriptorfile )
  
  return(distmat)
}



##check dist mat
.checkDistmat <- function( distmat, nodes, knodes = NULL ) {
  if( is.null(rownames(distmat)) | is.null(colnames(distmat)) ){
    stop("Error: no rownames or colnames of distmat.")
  }
  
  if( length(setdiff(nodes, rownames(distmat))) > 0 )
    stop("Error: some nodes not included in distmat")
  
  if( !is.null(knodes) ) {
    if( length(setdiff(knodes, rownames(distmat))) > 0 )
      stop("Error: some known nodes not included in distmat" )
  }
}






######################################
##get closeness
.DNA.closeness <- function( subDistMat, cpus = 1 ) {
  
  res <- apply(subDistMat, 1, sum, na.rm = TRUE)
  res <- 1/res
  res[which(res == Inf)] <- NA
  res[which(res == -Inf)] <- NA
  res
}





.DNA.eccentricity <- function( subDistMat, cpus = 1 ) {
  
  res <- apply(subDistMat, 1, max, na.rm = TRUE)
  res[which(res == Inf)] <- NA
  res[which(res == -Inf)] <- NA
  res
}




###################################################################################
##subgraph: distance 2 known stress genes
##from: index of nodes in the network
##to: index of known nodes in the network
####################################################################################
.DNA.distance2knowngenes <- function( subDistMat, mean = TRUE, cpus = 1 ) { 
  
  res <- NULL
  if( mean )
    res <- apply( subDistMat, 1, mean, na.rm = TRUE )
  else
    res <- apply( subDistMat, 1, sum, na.rm = TRUE)
  res 
}




#reference:Xu JZ and Li YJ. Bioinforamtics, 2006, 22(22): 2800-2805
.DNA.NIndex2knodes <- function( adjmat, connectivitylist, nodes, knodes, types = c("1N", "2N"), cpus = 1 ) {
  
  if( is.null(knodes) ) {
    stop("Error: knodes is not defined.")
  }
  
  refNodes <- rownames(adjmat)
  knodes.idx <- .getIndex( refNodes, knodes  )
  
  .fun_Index <- function( curNode, refNodes,  descfile, connectivitylist, knodes.idx, types ) {
    
    res <- rep(NA, 2)
    conn.idx <- connectivitylist[[curNode]]$all
    conn.len <- length(conn.idx)
    if( conn.len == 0 ){
      return(res)
    }
    res[1] <- length( intersect(conn.idx, knodes.idx))/conn.len  #1NIndex
    
    if( (length(types) == 1) & (types[1] == "1N") ){
      return (res)
    }
    
    sum2knode <- 0
    sum <- 0
    for( i in seq(length(conn.idx)) ){
      curnode.idx <- conn.idx[i]
      curRefNode <- refNodes[curnode.idx]
      curconn.idx <- connectivitylist[[curRefNode]]$all
      sum <- sum + length(curconn.idx)
      sum2knode <- sum2knode + length( intersect(curconn.idx, knodes.idx) )
    }
    res[2] <- sum2knode/sum
    
    res
  }#end fun
  
  if( cpus == 1 ) {
    result <- apply( matrix(nodes, ncol=1), 1, .fun_Index, refNodes = refNodes,connectivitylist = connectivitylist, knodes.idx = knodes.idx, types = types) 
  }else {
    sfInit(parallel = TRUE, cpus = cpus)
    sfLibrary("bigmemory", character.only=TRUE )
    result <- sfApply( matrix(nodes, ncol=1), 1, .fun_Index, refNodes = refNodes, connectivitylist = connectivitylist, knodes.idx = knodes.idx, types = types)
    sfStop()    
  }
  
  result <- t(result)
  rownames(result) <- nodes
  colnames(result) <- c("1N", "2N")
  result <- result[, types]  
  if( is.vector(result) )
    result <- matrix( result, ncol = 1)
  colnames(result) <- types
  
  result
}










