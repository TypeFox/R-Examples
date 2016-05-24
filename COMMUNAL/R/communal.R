## COMMUNAL package 1.1
## part 1: run clustering algorithms - see clusterRange() for a convenience  harness

makeCommunal <- setRefClass( "COMMUNAL", 
             fields = list( 
               cluster.list = "list",
               measures     = "array",
               clus.methods = "character",
               ks           = "numeric",
               validation   = "character",
               dist.metric  = "character", 
               item.names   = "character",
               call         = "call"
             ), 
             methods = list( 
               getClustering = function(k) {
                 # Returns the cluster assignment of each clustering method
                 #   for the given value of k.
                 # Args:
                 #   k: The desired number of clusters
                 #
                 # Returns:
                 #   Data frame of the cluster assignments, one column for each method.
                 if (!k %in% ks) {
                   stop(paste("k must be one of the original ks: ", 
                              paste(ks, collapse=", ")))
                 }
                 clusters <- eval(parse(text=paste("cluster.list$'", k, "'", sep="")))
                 df <- as.data.frame(clusters)
                 data.frame(lapply(df, as.integer), row.names = item.names)
               },
               show = function() {
                 cat("Reference class", classLabel(class(.self)), "\n\n")
                 cat("\tks:", paste(ks, collapse=", "), "\n\n")
                 cat("\tCall: ")
                 methods::show(call)
                 cat("\n(Use the getClustering method to extract data frame of cluster assignments)")
               }
               )
             )




## version 1.1 Added gap statistic functionality, and parallel options
## also verbosity
COMMUNAL <- function (data, ks, clus.methods = c("hierarchical", "kmeans", "diana",
                                                 "som", "sota", "pam", "clara", "agnes"), 
                      validation = c("Connectivity", "dunn", "wb.ratio", "g3", 
                                     "g2", "pearsongamma", "avg.silwidth", "sindex"), 
                      dist.metric = "euclidean", aggl.method = "ward", 
                      neighb.size = 10, seed = NULL, parallel=F, gapBoot=20, 
                      verbose=F, mc.cores=NULL, ...) 
{
  data <- as.matrix(data)
  if (is.null(colnames(data))) {
    colnames(data) <- paste("Sample", 1:ncol(data), sep = "")
    cat(paste("Added column names to data, from Sample1 to Sample", 
              ncol(data), "\n", sep = ""))
  }
  clus.methods = unique(match.arg(clus.methods, c("hierarchical", "kmeans", "diana",
                                                  "fanny", "som", "model", "sota", 
                                                  "pam", "clara", "agnes", "ccp-hc",
                                                  "ccp-km", "ccp-pam", "nmf"),
                                  several.ok = TRUE))
  
  possMeasures <- c("Connectivity", "average.between",
                    "g2", "ch", "sindex", "avg.silwidth", 
                    "average.within", "dunn", "widestgap",
                    "wb.ratio", "entropy", "dunn2", 
                    "pearsongamma", "g3", "within.cluster.ss", 
                    "min.separation", "max.diameter", 
                    "gapStatistic")
  if("all" %in% validation){
    validation <- possMeasures
  } else {
    validation <- unique(match.arg(validation, possMeasures, several.ok = TRUE))
  }
  
  if(("model" %in% clus.methods) & prod(dim(data)) > 50000){
      warning("Model-based clustering with mclust causes memory issues with large objects")
  }
  
  dist.metric = match.arg(dist.metric, c("euclidean", "correlation", 
                                         "manhattan"))
  aggl.method = match.arg(aggl.method, c("ward", "ward.D", 
                                         "ward.D2", "single", "complete", "average"))
  cl.methods = intersect(clus.methods, c("hierarchical", "kmeans", "diana", 
                                         "fanny", "som", "model", "sota", 
                                         "pam", "clara", "agnes"))
  ccp.methods = intersect(clus.methods, c("ccp-hc", "ccp-km", 
                                          "ccp-pam"))
  if ("nmf" %in% clus.methods) {
    if (!requireNamespace("NMF", quietly = TRUE)) {
      stop("NMF not installed.")
    }
    if (min(data) < 0) {
      stop("NMF requires all values in 'data' to be positive.")
    }
  }
  if (length(ccp.methods) > 0 && !requireNamespace("ConsensusClusterPlus", 
                                                   quietly = TRUE)) {
    stop(paste("Package ConsensusClusterPlus is required for clustering with", 
               paste(ccp.methods, collapse = ", ")))
  }
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
    abs(x - round(x)) < tol
  }
  if (!all(sapply(ks, function(x) is.wholenumber(x) & x > 
                    1))) {
    stop(paste("'ks' is invalid. Must be whole numbers > 1:", 
               paste(ks, collapse = ", ")))
  }
  if (any(duplicated(ks))) {
    stop("ks is invalid. Must be unique.")
  }
  getClus <- function(m, obj, k) {
    switch(m, hierarchical = cutree(obj, k), 
           diana = cutree(obj, k), 
           agnes = cutree(obj, k), 
           kmeans = obj[[as.character(k)]][["cluster"]], 
           clara = obj[[as.character(k)]][["clustering"]], 
           pam = obj[[as.character(k)]][["clustering"]], 
           fanny = obj[[as.character(k)]][["clustering"]], 
           sota = obj[[as.character(k)]][["clust"]], 
           model = obj[[as.character(k)]][["classification"]], 
           som = obj[[as.character(k)]][["unit.classif"]])
  }
  cluster.list <- replicate(length(ks), list())
  names(cluster.list) <- as.character(ks)
  mes <- array(dim = c(length(validation), length(ks), 0), 
               dimnames = list(c(), ks, c()))
  extractArgs <- function(params, valid.params) {
    if (is.list(valid.params)) {
      valid.params <- names(valid.params)
    }
    params[names(params) %in% valid.params]
  }
  t.data <- t(data)
  if (dist.metric == "correlation") {
    distance <- as.dist(1 - cor(data, use = "pairwise.complete.obs"))
  }
  else {
    if(verbose) cat("\nCalculating distance matrix... \n")
    distance <- dist(t.data, method = dist.metric)
  }
  if (length(cl.methods) > 0) {
    addl <- extractArgs(list(...), formals(clValid::clValid))
    if(verbose) cat("Clustering ", unlist(cl.methods), "... \n")
    if (length(addl) > 0) {
      cl <- clValid(t.data, ks, clMethods = cl.methods, 
                    validation = "internal", metric = dist.metric, 
                    method = aggl.method, neighbSize = neighb.size, 
                    verbose=verbose, addl)
    } else {
      cl <- clValid(t.data, ks, clMethods = cl.methods, 
                    validation = "internal", metric = dist.metric, 
                    method = aggl.method, neighbSize = neighb.size, 
                    verbose=verbose)
    }
    cobjs <- cl@clusterObjs
    for (k in ks) {
      k.clusts <- replicate(length(cl.methods), list())
      names(k.clusts) <- cl.methods
      for (m in cl.methods) {
        k.clusts[[m]] <- getClus(m, cobjs[[m]], k)
      }
      cluster.list[[as.character(k)]] <- k.clusts
    }
    for (alg in cl.methods) {
      if(verbose) cat("Calculating validation metrics (\".\" for each K) ")
      layer <- lapply(ks, function(k){
          if(verbose) cat(".")
          clus <- cluster.list[[as.character(k)]][[alg]]
          getValidation(clusters = clus,  data=data,
                        distance = distance, metric = dist.metric, 
                        validation = validation, neighb.size = neighb.size, 
                        num.clusters = k, alg = alg)      
      })
      layer <- Reduce(cbind, layer)
      colnames(layer) <- ks
      if(verbose) cat(" done \n")
      mes <- addLayer(mes, layer, name = alg)
    }
  }
  ## add Mclust?
  if (length(ccp.methods) > 0) {
    if (dist.metric == "manhattan") {
      warning("For ConsensusClusterPlus, 'manhattan' distance metric not available.\n
              Switching to 'euclidean' instead.")
    }
    d.metric <- switch(dist.metric, correlation = "pearson", 
                       manhattan = "euclidean", euclidean = "euclidean")
    addl <- extractArgs(list(...), formals(ConsensusClusterPlus::ConsensusClusterPlus))
    if (!"reps" %in% names(addl)) {
      addl$reps <- 500
    }
    for (alg in ccp.methods) {
      if(verbose) cat("Calculating ", alg, "  \n")
      alg1 <- sub("ccp-", "", alg)
      layer <- matrix(NA, nrow = dim(mes)[1], ncol = dim(mes)[2])
      rcc <- do.call(ConsensusClusterPlus::ConsensusClusterPlus, 
                     c(list(data, maxK = max(ks), clusterAlg = alg1, 
                            distance = d.metric, seed = seed), addl))
      if(verbose) cat("Calculating validation metrics (\".\" for each K) ")
      for (i in 1:length(ks)) {
        k = ks[i]
        clus <- rcc[[k]][["consensusClass"]]
        cluster.list[[as.character(k)]][[alg]] <- clus
      }
      layer <- lapply(ks, function(k){
          if(verbose) cat(".")
          clus <- rcc[[k]][["consensusClass"]]
          getValidation(clusters = clus,  data=data,
                        distance = distance, metric = dist.metric, 
                        validation = validation, neighb.size = neighb.size, 
                        num.clusters = k, alg = alg)      
      })
      layer <- Reduce(cbind, layer)
      colnames(layer) <- ks
      if(verbose) cat(" done \n")
      mes <- addLayer(mes, layer, name = alg)
    }
  }
  if ("nmf" %in% clus.methods) {
    layer <- list()
    addl <- extractArgs(list(...), formals(NMF::nmf))
    if(verbose) cat("Calculating NMF and validation metrics (\".\" for each K) ")
    for (i in 1:length(ks)) {
      k = ks[i]
      if (length(addl) > 0) {
        nmfRes <- NMF::nmf(data, rank = k, seed = seed, 
                           .options = "v", addl)
      }
      else {
        nmfRes <- NMF::nmf(data, rank = k, seed = seed, 
                           .options = "v")
      }
      clus <- as.numeric(NMF::predict(nmfRes))
      cluster.list[[as.character(k)]][["nmf"]] <- clus
      layer[[i]] <- getValidation(clusters = clus,  data=data,distance = distance, 
                                  metric = dist.metric, validation = validation, 
                                  neighb.size = neighb.size, num.clusters = k, 
                                  alg = "NMF")
      if(verbose) cat(".")
    }
    layer <- Reduce(cbind, layer)
    colnames(layer) <- ks
    if(verbose) cat(" done \n")
    mes <- addLayer(mes, layer, name = "nmf")
  }
  
  ## prevents runaway recursion in gap statistic calculation
  calls <- unlist(lapply(sys.calls(), function(x) 
                    strsplit(as.character(x), split="\\(")[[1]][1]))
  # cat(calls)
  if("gapStat_BootData" %in% calls){
      return(mes)
  }
  if("gapStatistic" %in% validation){
    gapOut <- gapStat_BootData(data, ks=ks, gapBoot=gapBoot,  
                               clus.methods = clus.methods, dist.metric = dist.metric, 
                               aggl.method = aggl.method, neighb.size = neighb.size,
                               seed = seed, parallel=parallel, mc.cores=mc.cores, ...)

    ## mod to subtract these from mes object
    for(alg in dimnames(mes)[[3]]){
      E.logW <- gapOut[, alg]
      mes["gapStatistic", , alg] <- E.logW - mes["gapStatistic", , alg]
    }
  }
  
  makeCommunal$new(cluster.list = cluster.list, measures = mes, 
                   clus.methods = clus.methods, ks = ks, validation = validation, 
                   dist.metric = dist.metric, item.names = colnames(data), 
                   call = match.call())
}


## v1.1: calculate 
gapStat_BootData <- function (data, ks, gapBoot=20, gapVerbose=F, 
                              clus.methods = c("hierarchical", "kmeans"), 
                              dist.metric = "euclidean", aggl.method = "average", 
                              neighb.size = 10, seed = NULL, parallel=F, 
                              mc.cores=NULL, ...) { 
  t.data <- t(data)
  
  n <- nrow(t.data)
  xs <- scale(t.data, center = TRUE, scale = FALSE)
  m.x <- rep(attr(xs, "scaled:center"), each = n)
  V.sx <- svd(xs)$v
  rng.x1 <- apply(xs %*% V.sx, 2, range)
  
  if(!is.null(seed)) set.seed(seed)
  
  cat("Gap statistic bootstrapping, ", gapBoot,
        " reps,  one \".\" per rep:\n")
  cat("(Each bootstrap recursively runs COMMUNAL)\n")
  
  scalePos <- F
  if("nmf" %in% clus.methods){
    scalePos <- T
    cat("Scaling all gap bootstraps to positive values to comply with NMF\n")
  }
  
  make_z <- function(rng.x1, n, V.sx, m.x, data, colnamesData, scalePos){
    z1 <- apply(rng.x1, 2, function(M) runif(n, min = M[1], max = M[2]))
    z <- tcrossprod(z1, V.sx) + m.x
    z <- t(z)
    colnames(z) <- colnamesData
    if(scalePos) z <- z - min(z) + 1
    z
  }
  
  logWks <- array(dim = c(gapBoot, length(ks), length(clus.methods)), 
                  dimnames = list(paste("gapBoot", seq_len(gapBoot)), 
                                  ks, clus.methods))
  colnames(logWks) <- paste0("K_", ks)
  
  if(parallel){
    ### will FAIL on Windows!!!!
    if(is.null(mc.cores)) mc.cores <- parallel::detectCores()-1
    out <- parallel::mclapply(1:gapBoot, mc.cores=mc.cores, mc.preschedule = FALSE,
                    mc.cleanup = FALSE, function(b){
        set.seed(b)
        z <- make_z(rng.x1, n, V.sx, m.x, data, colnames(data), scalePos)
         
        ## run each bootstrapped data object through COMMUNAL again, 
        ## but only eval gap statistic
        cat(".", if (b%%5 == 0) paste(b, "\n"))
        COMMUNAL(data=z, validation="gapStatistic", ks=ks, 
                    clus.methods=clus.methods, dist.metric=dist.metric, 
                    aggl.method=aggl.method, neighb.size = neighb.size, 
                    verbose=gapVerbose, ...)      
    })
    #reassign output into preallocated mtx
    for(b in 1:gapBoot){
        logWks[b, , ] <- out[[b]]
    }
  } else {
    for (b in 1:gapBoot) {
      set.seed(b)
      z <- make_z(rng.x1, n, V.sx, m.x, data, colnames(data), scalePos)
      
      ## run each bootstrapped data object through COMMUNAL again, 
      ## but only eval gap statistic
      logWks[b, , ] <- COMMUNAL(data=z, validation="gapStatistic", ks=ks, 
                                clus.methods=clus.methods, dist.metric, 
                                aggl.method, neighb.size = neighb.size, 
                                seed = seed, verbose=gapVerbose, ...)
      
      cat(".", if (b%%5 == 0) paste(b, "\n"))
    }
  }
  ##applying colmeans over MAR=3 returns algs in cols, K in rows
  E.logW <- apply(logWks, 3, colMeans)
  dim(E.logW) <- dim(logWks)[2:3] ## in case only one k, drops to vector
  dimnames(E.logW) <- dimnames(logWks)[2:3]
  # SE.sim <- sqrt((1 + 1/gapBoot) * apply(logWks, ?, function(x) apply(x, 2, var)))
  # gap = E.logW - logW
  return(E.logW)
}


## version 1.1 added gapStatistic functionality
getValidation <- function(clusters, alg, data = NULL, distance = NULL, 
                          metric = "euclidean", neighb.size = 10, 
                          validation = c("dunn", "avg.silwidth", "Connectivity"), 
                          num.clusters = -1) {
    # Computes the validation scores for a given clustering of the data.
    # Args:
    #   clusters: The cluster assignments
    #   alg: the algorithm name, used to print error messages
    #   data: Data matrix. Must specify either this or distance
    #   distance: Distance matrix
    #   metric: Used to calculate distance matrix
    #   neighb.size: used to calculate connectivity
    #   validation: types of validation metrics to evaluate
    #   num.clusters: the intended number of clusters
    # Returns:
    #   Data frame of the cluster assignments, one column for each method.
    result <- numeric()
    metric <- match.arg(metric, c("euclidean", "correlation", "manhattan"))
    if (all(is.na(clusters))) {
        warning(paste0(alg, ": no cluster assignments available - 
                       validation measures will be NA for k=", num.clusters))
        return(rep(NA, length(validation)))
    }
    if("gapStatistic" %in% validation){
        ##need to modify getValidation calls to pass data
        ii <- seq_len(ncol(data))
        
        ## calculate W for a given k
        W.k <-  sum(vapply(split(ii, clusters), function(I) {
            xs <- data[ , I, drop = FALSE]
            sum(dist(t(xs))/ncol(xs))
        }, 0))
        W.k <- log(0.5*W.k)
        ## if only gap, skip the rest
        if(length(validation)==1)  return(t(t(c("gapStatistic"=W.k))))
        result <- c(result, "gapStatistic"=W.k)
    }
    if (is.null(distance)) {
        if (is.null(data)) {
            stop("Missing data and distance args. Must specify at least one.")
        }
        distance <- dist(t(data), method=metric)
    }
    if ("Connectivity" %in% validation) {
        result <- c(result, "Connectivity" = clValid::connectivity(distance = distance, 
                                                                   clusters = clusters, 
                                                                   neighbSize = neighb.size))
    }
    run.g2 <- "g2" %in% validation
    run.g3 <- "g3" %in% validation
    run.sepindex <- "sindex" %in% validation
    ## can't compute silhouette if only one cluster
    ifelse(("avg.silwidth" %in% validation) & (length(unique(clusters))>1), 
           run.silhouette <- T,
           run.silhouette <- F)
    stats <- unlist(fpc::cluster.stats(d = distance, clustering = clusters, 
                                       silhouette = run.silhouette, G2 = run.g2, 
                                       G3 = run.g3, sepindex = run.sepindex))
    
    if("avg.silwidth" %in% validation){
        validation <- validation[validation != "avg.silwidth"]
        result <- c(result, stats[validation[validation %in% names(stats)]])
        ifelse(run.silhouette, 
               result <- c(result, stats["avg.silwidth"]),
               result <- c(result, "avg.silwidth" = NA))
    } else {
        result <- c(result, stats[validation[validation %in% names(stats)]])
    }
    t(t(result))
}

addLayer <- function(arr, layer, name) {
    for(i in 1:2){
        dimnames(arr)[[i]] <- dimnames(layer)[[i]] 
    }    
    # Adds another layer in the third dimension of the array.    
    new.dim <- dim(arr)
    new.dim[3] <- new.dim[3] + 1
    new.names <- dimnames(arr)
    if (dim(arr)[3] == 0) {
        new.names[[3]] <- name
    } else {
        new.names[[3]] <- c(new.names[[3]], name)
    }
    new.array <- array(c(arr,as.vector(layer)), dim = new.dim, dimnames = new.names)
    if (dim(arr)[3] > 0) {
        for (i in 1:dim(arr)[3]) {
            stopifnot(identical(arr[, , i], new.array[, , i]))
            
        }
    }
    new.array
}
