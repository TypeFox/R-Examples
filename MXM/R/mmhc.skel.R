mmhc.skel <- function(dataset, max_k = 3, threshold = 0.05, test = NULL, rob = FALSE, fast = FALSE, nc = 1, graph = FALSE) {
  ## dataset is either conitnuous or categorical data  
  ## max_k is the maximum number of variables upon which to condition
  ## threshold is the level of significance to reject the independence
  ## test can be either testIndFisher (default) or testIndSpearman for continuous data
  ## OR gSquare (default) for categorical data
  ## rob is for robust correlation
  ## nc is the number of cores to use, set to 1 by default
  
  dataset <- as.matrix(dataset)
  n <- ncol(dataset)
  G <- matrix(0, n, n )
  if (test == "testIndSpearman") {
    dataset <- apply(dataset, 2, rank)
    rob = FALSE 
  }
  
  if (nc == 1 || nc < 1 || is.null(nc) ) {
    if (fast == TRUE) {   
      pa <- proc.time()
      a <- MMPC(1, dataset, max_k = max_k, test = test, threshold = threshold, robust = rob)
      sel <- a@selectedVars
      G[1, sel] <- 1 
    
      for ( i in 2:n ) {
        ina <- 1:n
        che <- which( G[ 1:c(i - 1), i ] == 0 ) 
        ina[c(i, che)] <- 0   ;  ina <- ina[ina>0]
        a <- MMPC(dataset[, i], as.matrix(dataset[, -c(i, che)]), max_k = max_k, test = test, threshold = threshold, robust = rob)
        if ( !is.null(a) ) {
          sel <- a@selectedVars
          sela <- ina[sel]
          G[i, sela] <- 1 
        }  else {
           G[i, ] <- 0
        }
      }
      runtime <- proc.time() - pa
      
    } else {
      pa <- proc.time()
      for (i in 1:n) {
        a <- MMPC(i, dataset, max_k = max_k, test = test, threshold = threshold, robust = rob)
        sel <- a@selectedVars
        G[i, sel] <- 1 
      } 
      runtime <- proc.time() - pa
    }
  }  else {
   
    pa <- proc.time() 
    cl <- makePSOCKcluster(nc)
    registerDoParallel(cl)
    sel <- numeric(n)
    mod <- foreach(i = 1:n, .combine = rbind, .export = c("MMPC") ) %dopar% {
      ## arguments order for any CI test are fixed
      sel <- numeric(n)
      a <- MMPC(i, dataset, max_k = max_k, test = test, threshold = threshold, robust = rob)
      sel[a@selectedVars] <- 1
      return(sel)
    }
    stopCluster(cl)
    G <- as.matrix(mod)
    runtime <- proc.time() - pa
  }
   
  G2 <- G - t(G)
  G[ G2 != 0 ] <- 0
  diag(G) <- 0
  
  info <- summary( rowSums(G) )
  density <- sum(G) / ( n * ( n - 1 ) )
  
  if (is.null( colnames(dataset) ) ) {
    colnames(G) <- rownames(G) <- paste("X", 1:n, sep = "")
  } else  colnames(G) <- rownames(G) <- colnames(dataset)
  
  if(graph == TRUE)
  {
    if(requireNamespace("Rgraphviz", quietly = TRUE, warn.conflicts = FALSE) == TRUE)
    {
      am.graph <- new("graphAM", adjMat = G, edgemode = "undirected")
      plot( am.graph, main = paste("Skeleton of the MMHC algorithm for", deparse( substitute(dataset) ) ) )
    }else{
      warning('In order to plot the generated network, package Rgraphviz is required.')
    }
  }
  
  list(runtime = runtime, density = density, info = info, G = G)
}