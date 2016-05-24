######################################################################################
##from gene expression to transcriptional network
exp2net <- function( expmat, method = c("GCC", "PCC", "SCC", "KCC", "BiWt", "MI", "MINE"),
                     pvalue = 0.01, cpus = 1, expDescribe = "Control", connListFlag = TRUE, distmatFlag = TRUE, saveType = "bigmatrix", netResFileDic, ...  ) {

     call <- match.call()

    if( is.null(rownames(expmat) ) )
      stop("Error: no rownames for expmat.\n")

    if( length(method) > 1 )
      method <- method[1]

    ##create PSOLResDic
    dir.create( path = netResFileDic, showWarnings = FALSE)

    #if( !require(rsgcc) ) {
    #  install.packages("rsgcc")
    #  library(rsgcc)
    #}
 
    ##get adjmatrix 
    cat("...calculating adjacency matrix...\n")
    backingpath = netResFileDic
    backingfile <- paste( expDescribe, "_", method, "_adjmat_bfile", sep = "" )
    descriptorfile <- paste( expDescribe, "_", method, "_adjmat_dfile", sep = "" )
    adjmat <- adjacencymatrix( mat = expmat,  method = method, cpus = cpus, saveType= saveType, 
                              backingpath= backingpath, backingfile= backingfile, descriptorfile= descriptorfile )

   
    ##get threshold
    cat("...calculating correlation coefficient at the given significance level...\n")
    backingfile_padjmat <- paste( expDescribe, "_", method, "_permut_adjmat_bfile", sep = "" )
    descriptorfile_padjmat <- paste( expDescribe, "_", method, "_permut_adjmat_dfile", sep = "" )
    threshold <- .cor.threshold( expMat = expmat, sigLevels = pvalue, corMethod = method, distMethod = "permutation", tailed = "two.sided", cpus = cpus,
                                 backingpath= backingpath, backingfile= backingfile_padjmat, descriptorfile= descriptorfile_padjmat  )

    ##get igraph project
    cat("...build network...")
    mygraph <- .formatAdj2Graph( adjmat = adjmat, threshold = threshold, file = paste(backingpath, expDescribe, "_", method, "_graph", sep = "" ), format = c("edgelist" ) )

    ##connectionList
    connectivityList <- NULL
    if( connListFlag == TRUE ){
       cat( "...calculating connectivity List...\n")
       connectivityList <- .DNA.connectivitylist( adjmat = adjmat, threshold = threshold, backingpath = netResFileDic, descriptorfile = descriptorfile, nodes = rownames(expmat), cpus = cpus )  
       save( connectivityList, file = paste( netResFileDic, expDescribe, "_net_connectivityList.RData", sep = "") )
    }

    ##distance matrix
    distmatrix <- NULL
    if( distmatFlag == TRUE ) {
        if( !exists("v") ){
          v <- V(mygraph)$name
        }else {
          v <- intersect(v, V(mygraph)$name )
        }
        if( !exists("to") ){
          to <- V(mygraph)$name
        }else {
          to <- intersect(to, V(mygraph)$name )
        }
       
       cat("...calculating distance matrix...\n")
       distmatrix <- .DNA.distance( graph = mygraph, v = v, to = to, cpus = cpus, saveType = saveType, 
                                    backingpath = netResFileDic, 
                                    backingfile = paste( expDescribe, "_", method, "_distmat_bfile", sep = "" ), 
                                    descriptorfile = paste( expDescribe, "_", method, "_distmat_dfile", sep = "" ) )
  
    }

    res <- list( expmat = expmat, method = method, pvalue = pvalue, expDescribe = expDescribe, netResFileDic = netResFileDic,
                 adjmat = adjmat, adjmat_backingfile = backingfile, adjmat_descriptorfile = descriptorfile, 
                 threshold = threshold, graph = mygraph, connectivityList = connectivityList, distmatrix = distmatrix )
    save( res, file = paste( netResFileDic, expDescribe, "_network_exp2net.RData", sep = "" ) )
    res
}






###generate network feature matrix with user-specified network characteristics.
netFeatureMatrix <- function( net1, net2, nodes = NULL, knodes = NULL, cpus = 1, verbose = TRUE, netResFileDic, 
                              features = c( "expDistance", "ASC", "corDistance", "AllConnectivity", 
                                            "PosConnectivity", "NegConnectivity", "closeness", "eccentricity", 
                                            "eigenvector", "page.rank", "dis2knodes", "closeness2knodes", "eccenticity2knodes") ) {


   ##check features
   if( length(features) == 0 ) 
      stop("Error: no features is specified.\n" )

   tmp <- setdiff( features, c( "expDistance", "ASC", "corDistance", "AllConnectivity", "PosConnectivity", "NegConnectivity", "closeness", "eccentricity", "eigenvector", "page.rank", "dis2knodes", "closeness2knodes", "eccenticity2knodes") )
   if( length(tmp) > 0 ) {
     cat(tmp)
     stop("\nError: these features not defined in the current version of mlDNA.\n")
   }
   ##nodes in network
   graphNodes1 <- V(net1$graph)$name
   graphNodes2 <- V(net2$graph)$name
   interSectNodes <- intersect( graphNodes1, graphNodes2 )
  
  ##if nodes not specified, all nodes will be considered
  if( is.null(nodes) ) {
    nodes <- interSectNodes
  }else {
    nodes <- intersect( interSectNodes, nodes )
  }
  
  ##get connectivityList
  connectivityList1 <- net1$connectivityList
  if( is.null(connectivityList1) ) {
     if(verbose)
       cat("...calculating connectivityList1...\n")
     connectivityList1 <- .DNA.connectivitylist( adjmat = net1$adjmat, threshold = net1$threshold, backingpath = net1$netResFileDic, descriptorfile = net1$adjmat_descriptorfile, nodes = nodes, cpus = cpus )  
     save( connectivityList1, file = paste( net1$netResFileDic, net1$expDescribe, "_net_connectivityList.RData", sep = "") )
  }
  connectivityList2 <- net2$connectivityList
  if( is.null(connectivityList2) ) {
      if(verbose)
       cat("...calculating connectivityList2...\n")
     connectivityList2 <- .DNA.connectivitylist( adjmat = net2$adjmat, threshold = net2$threshold, backingpath = net2$netResFileDic, descriptorfile = net2$adjmat_descriptorfile, nodes = nodes, cpus = cpus )  
     save( connectivityList2, file = paste( net2$netResFileDic, net2$expDescribe, "_net_connectivityList.RData", sep = "") )
  }

  ##start to clacluate properties
  leftFeatures <- features
  propmat <- NULL

  ##for three network difference.
  difFeatures <- intersect( leftFeatures, c("expDistance", "corDistance", "ASC" ) )
  if( length(difFeatures) > 0 ) {
     matd <- .network.differences( expmat1 = net1$expmat, net1 = net1$adjmat, threshold1 = net1$threshold, backingpath1 = net1$netResFileDic, descriptorfile1 = net1$adjmat_descriptorfile, connectivityList1 = connectivityList1,  
                                   expmat2 = net2$expmat, net2 = net2$adjmat, threshold2 = net2$threshold, backingpath2 = net2$netResFileDic, descriptorfile2 = net2$adjmat_descriptorfile, connectivityList2 = connectivityList2,  
                                   nodes = nodes, properties = difFeatures, cpus = cpus, verbose = verbose ) 
    colnames(matd) <- difFeatures
    leftFeatures <- setdiff( leftFeatures, difFeatures )
  }
  if( length(leftFeatures) == 0 ) {
     return(matd)
  }


  ##further network centralities
  mat1 <- .network.properties(nodes = nodes, adjmat = net1$adjmat, threshold = net1$threshold,  backingpath = net1$netResFileDic, backingfile = net1$adjmat_backingfile, descriptorfile = net1$adjmat_descriptorfile,  
                               graph = net1$graph, distmat = net1$distmatrix, connectivityList = connectivityList1, knodes = knodes, cpus = cpus, 
                               properties = leftFeatures, netDescribe = net1$expDescribe, verbose = verbose )
  
  mat2 <- .network.properties(nodes = nodes, adjmat = net2$adjmat, threshold = net2$threshold,  backingpath = net2$netResFileDic, backingfile = net2$adjmat_backingfile, descriptorfile = net2$adjmat_descriptorfile,  
                               graph = net2$graph, distmat = net2$distmatrix, connectivityList = connectivityList2, knodes = knodes, cpus = cpus, 
                               properties = leftFeatures, netDescribe = net2$expDescribe, verbose = verbose )
  mat2 <- mat2[, colnames(mat1)]

  if( (nrow(mat1) != nrow(mat2) ) | (ncol(mat1) != ncol(mat2) ) ) {
     stop("Error: different dimensions for feature matrix from net1 and net2.\n" )
  }
  matdif <- mat1 - mat2[rownames(mat1), ]
  colnames(matdif) <- paste( colnames(matdif), ".d", sep = "" )
  colnames(mat1) <- paste( colnames(mat1), net1$expDescribe, sep = "." )
  colnames(mat2) <- paste( colnames(mat2), net2$expDescribe, sep = "." )

  resmat <- cbind( matd, mat1, mat2, matdif )

  resmat

}

interactionRemoval <- function( adjmat1, adjmat2, threshold1, threshold2 ) {

  getSetDiff <- function( v1, v2 ) {
     len1 <- length(v1)
     if( len1 < 1e+09 ){
       return( setdiff(v1, v2) )
     }else {
       idx1 <- floor( len1/2 )
       res1 <- setdiff( v1[1:idx1], v2)
       res2 <- setdiff( v1[(idx1+1):len1], v2)
       return( c(res1, res2))
    }
  }#end fun

  #check dim
  if( nrow(adjmat1) != nrow(adjmat2) ) {
     stop("Error: different genes in adjmat1 and adjmat2!")
  }

  if( rownames(adjmat1) == NULL | rownames(adjmat2) == NULL ) {
     stop("Error: no rownames assigned for adjmat1 or adjmat2!")
  }
 
  #check gene order
  if(  all.equal( rownames(adjmat1), rownames(adjmat2) ) != TRUE ) {
     stop("Error: different order of genes in adjmat1 and adjmat2!")
  }

  #now start
  mat <- matrix(0, nrow= nrow(adjmat1), ncol=1 )
  rownames(mat) <- rownames(adjmat1)
  colnames(mat) <- c("net2_Specific_Degree")

 
  net1_Idx <- which( abs(adjmat1[,]) >= threshold1 )
  net2_Idx <- which( abs(adjmat2[,]) >= threshold2 )
  idx <- getSetDiff( net2_Idx, net1_Idx)
  
  if( length(idx) != 0 ){
    tmp <- adjmat2[,]
    tmp[idx] <- 1
    netS_Idx_0 <- which( abs(adjmat2[,]) < threshold2 )
    tmp[netS_Idx_0] <- 0
    tmp[net1_Idx] <- 0
    mat[,1] <- apply( tmp, 1, sum )
  }
  mat
}



