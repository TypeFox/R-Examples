clusterTree <-
function(X, k, h = NULL, density = "knn", dist = "euclidean", d = NULL,
         Nlambda = 100, printProgress = FALSE) {

  if (!is.numeric(X) && !is.data.frame(X)) {
    stop("X should be an n by d matrix of coordinates")
  }
  if (!is.numeric(k) || length(k) != 1 || k < 1) {
    stop("k should be a positive integer")
  }
  if (!is.null(h) && (!is.numeric(h) || length(h) != 1)) {
    stop("h should be a real value")
  }
  if (density != "knn" && density != "kde") {
    stop("density should be a string: either 'knn' or 'kde'")
  }
  if (dist != "euclidean" && dist != "arbitrary") {
    stop("dist should be a string: either 'euclidean' or 'arbitrary'")
  }
  if (density == "kde" && dist != "euclidean") {
    stop("kde is only possible with dist = 'euclidean'")
  }
  if (!is.null(Nlambda) && (!is.numeric(Nlambda) || length(Nlambda) != 1)) {
    stop("Nlambda should be a number")
  }
  if (!is.logical(printProgress)) {
    stop("printProgress should be logical")
  }
  if (!is.null(d) && (!is.numeric(d) || length(d) != 1)) {
    stop("d should be a number")
  }
      
  if (dist == "euclidean") { 
    n <- dim(X)[1]
    d <- dim(X)[2]
  } else if (dist == "arbitrary") { 
    n <- dim(X)[1]
    if (is.null(d)) {
      d <- 1
    }
    distMat <- X
  }


  ## start adjacency matrix
  adjMat <- diag(n)

  ## Compute density estimator: knn or kde
  if (density == "knn" && dist == "euclidean") {
    knnInfo <- FNN::get.knn(X, k = k, algorithm = "kd_tree")
    for(i in seq_len(n)) {
      adjMat[i, knnInfo[["nn.index"]][i, ]] <- 1
    }
    r.k <-  apply(knnInfo[["nn.dist"]], 1, max)   
    v.d <- pi ^ (d / 2) / gamma(d / 2 + 1)
    hat.f <- k / (n * v.d * r.k ^ d)  
  } else if (density == "knn" && dist == "arbitrary") {
    orderMat <- apply(distMat, 2, order)[2:(k+1), ]
    for(i in seq_len(n)) {
      adjMat[i, orderMat[, i]] <- 1   
    }
    r.k <- apply(apply(X, 2, sort)[2:(k+1), ], 2, max)
    v.d <- pi ^ (d / 2) / gamma(d / 2 + 1)
    hat.f <- k / (n * v.d * r.k ^ d)  
  } else if (density == "kde") {  #kde estimate
    knnInfo <- FNN::get.knn(X, k = k, algorithm = "kd_tree")
    for(i in seq_len(n)) {
      adjMat[i, knnInfo[["nn.index"]][i, ]] <- 1
    }
    hat.f <- kde(X, X, h) 
  }

  # ordered value of the density
  ord.hat.f <- order(hat.f)

  # starting graph  
  G <- igraph::graph.adjacency(adjMat, mode = "undirected")   ## could be changed to directed

  # Lambda grid
  Lambda <- hat.f[ord.hat.f]
  if (is.null(Nlambda)) { 
    Nlambda <- n
  } else {
    Nlambda <- min(n, Nlambda)
    Lambda <- seq(min(Lambda), max(Lambda), length = Nlambda)
  }

  exclude <- numeric()

  ## in CLUSTERS we store the clusters found for each level of lambda_j
  CLUSTERS <- list()
  if (printProgress) {
    cat("0   10   20   30   40   50   60   70   80   90   100\n")
    cat("|----|----|----|----|----|----|----|----|----|----|\n")
    cat("*")    
  }
  percentageFloor <- 0
  for (j in seq_len(Nlambda)) {
    OldExcluded <- exclude
    lambda <- Lambda[j]
    present <- which(hat.f >= lambda)  
    exclude <- setdiff(seq_len(n), present)  # points with density less than lambda
    NewExcluded <- setdiff(exclude,OldExcluded)  # the new excluded point  
    G[NewExcluded, present] <- FALSE  # remove edges of the new excluded point
    clust <- igraph::clusters(G)
    CLUSTERS[[j]] <- list("no" = clust[["no"]], "mem" = clust[["membership"]],
                          "present" = present, "exclude" = exclude)   

    if (printProgress &&
        floor((100 * j / Nlambda - percentageFloor) / 2) > 0) {
      for (aa in seq_len(floor((100 * j / Nlambda - percentageFloor) / 2))) {
        cat("*")
        percentageFloor <- percentageFloor + 2
      }
    }
  }
  if (printProgress) {
    cat("\n")
  }

  ## Now assign ID, Generation and Components to each new cluster
  id <- 0                  # id assigned to each new cluster
  components <- list()     # data points contained in each new cluster
  generation <- numeric()  #generation of each new cluster
  for (j in seq_len(Nlambda)) {
    presentMembership <-
        unique(CLUSTERS[[j]][["mem"]][CLUSTERS[[j]][["present"]]])
    for (i in presentMembership) {
      id <- id + 1
      components[[id]] <- which(CLUSTERS[[j]][["mem"]] == i)
      generation[id] <- j       
    }
  }

  # find the father of each cluster
  father <- numeric()
  startF <- which(generation == 2)[1]
  for (i in startF:length(components)) {
    for (j in which(generation == (generation[i]-1))) {
      if (setequal(intersect(components[[i]], components[[j]]),
                   components[[i]])) {
        father[i] <- j
        break
      }
    }
  }
  father[is.na(father)] <- 0  #for the roots set father = 0

  ## Convert the clusters into BRANCHES
  bb <- 0               # count number of branches
  branch <- numeric()   # map each cluster into a branch
  base <- numeric()     # x coordinate of the base of each branch
  top <- numeric()      # y-top of each branch
  bottom <- numeric()   # y-bottom of each branch
  compBranch <- list()  # data points corresponding to this branch
  silo <- list()        # x coordinates of each branch silo
  rank <- numeric()     # rank among brothers
  parent <- numeric()   # father of each branch
  sons <- list()        # sons of each branch 

  # if there is more than 1 root create a fake single root
  if (sum(generation == 1) > 1) {
    bb <- bb + 1
    silo[[bb]] <- c(0, 1)
    base[bb] <- 0.5
    compBranch[[bb]] <- seq_len(n)
    rank[bb] <- 1
    parent[bb] <- 0  # parent of roots is set to be 0
    top[bb] <- 0
    bottom[bb] <- 0  # y-bottom of roots is set to be 0
  }
  for (i in seq(along = father)) {

    # the first generation is treated separately    
    # multiple roots are children of the fake single root
    if (sum(generation == 1) > 1 & generation[i] == 1) { 
      Bros <- which(generation == generation[i])  
      bb <- bb+1
      branch[i] <- bb
      rank[bb] <- sum(generation[seq_len(i)] == generation[i] &
                      father[seq_len(i)] == father[i])  #same gen, same father.
      silo[[bb]] <- siloF(c(0,1), length(Bros), rank[bb]) 
      base[bb] <- sum(silo[[bb]]) / 2
      top[bb] <- min(hat.f[components[[i]]])
      compBranch[[bb]] <- components[[i]]   
      parent[bb] <- 1  # parent of roots is set to be 1
      bottom[bb] <- 0  # y-bottom of roots is set to be 0
      ## add this branch to the list of sons of its parent
      if (length(sons) < parent[bb]) {
        sons[[parent[bb]]] <- bb
      } else {
        sons[[parent[bb]]] <- c(sons[[parent[bb]]], bb)
      }   
    } else if (sum(generation == 1) == 1 & generation[i] == 1){ #is there is 1 root
      bb <- bb + 1
      branch[i] <- bb
      silo[[bb]] <- c(0, 1)
      base[bb] <- 0.5
      top[bb] <- min(hat.f[components[[i]]])
      compBranch[[bb]] <- components[[i]]   
      parent[bb] <- 0   # parent of roots is set to be 0
      bottom[bb] <- 0   # y-bottom of roots is set to be 0
    } else {  
      Bros <- which(generation == generation[i] & father == father[i])  
      ## if the cluster has brothers, then there is a split and new branches are created
      if (length(Bros) > 1) {
        bb <- bb + 1
        branch[i] <- bb
        parent[bb] <- branch[father[i]]
        rank[bb] <- sum(generation[seq_len(i)] == generation[i] &
                        father[seq_len(i)] == father[i])
        silo[[bb]] <- siloF(silo[[parent[bb]]], length(Bros), rank[bb]) 
        base[bb] <- sum(silo[[bb]]) / 2
        top[bb] <- min(hat.f[components[[i]]]) 
        bottom[bb] <- top[parent[bb]]
        compBranch[[bb]] <- components[[i]]
        ## add this branch to the list of sons of its parent
        if (length(sons) < parent[bb]) {
          sons[[parent[bb]]] <- bb
        } else {
          sons[[parent[bb]]] <- c(sons[[parent[bb]]], bb)
        }     
      }

      ## if the cluster does not have brothers, then no new branches are created
      ## and this cluster is assigned to an old branch  
      if (length(Bros) == 1){  
        for (j in which(generation == (generation[i] - 1))) {
          if (setequal(intersect(components[[i]], components[[j]]),
                       components[[i]]))
          belongTo <- branch[j]
        }
        top[belongTo] <- min(hat.f[components[[i]]]) #update top of the branch
        branch[i] <- belongTo
      }
    }
  } 


  #info for alpha Tree
  ID <- seq_len(bb)
  alphaBottom <- numeric(bb)
  alphaTop <- numeric(bb)
  for (i in seq_len(bb)) {
    alphaBottom[i] <- sum(hat.f > bottom[i]) / n
    alphaTop[i] <- sum(hat.f > top[i]) / n
  }

  ## info for the kappa tree
  kTree <- findKtree(bb, parent, sons, compBranch, n)
  kappaTop <- kTree[["kappaTop"]]
  kappaBottom <- kTree[["kappaBottom"]]


  ## r Tree
  if (density != "kde") {
    rTop <- (k / (n * v.d * top)) ^ (1 / d)
    rBottom <- (k / (n * v.d * bottom)) ^ (1 / d)

    for (i in seq_len(bb)) {
      if (bottom[i] == 0) {
        rBottom[i] <- max(r.k)
      }
      if (top[i] == 0) {
        rTop[i] <- max(r.k)
      }
    }

    out <- list("density" = hat.f, "DataPoints" = compBranch, "n" = n,
           "id" = seq_len(bb), "sons" = sons, "parent" = parent, "silo" = silo,
           "Xbase" = base, "lambdaBottom" = bottom, "lambdaTop" = top,
           "rBottom" = rBottom, "rTop" = rTop,
           "alphaBottom" = alphaBottom, "alphaTop" = alphaTop,
           "kappaBottom" = kappaBottom, "kappaTop" = kappaTop)
  } else {
    out <- list("density" = hat.f, "DataPoints" = compBranch, "n" = n,
           "id" = seq_len(bb), "sons" = sons, "parent" = parent, "silo" = silo,
           "Xbase" = base, "lambdaBottom" = bottom, "lambdaTop" = top,
           "alphaBottom" = alphaBottom, "alphaTop" = alphaTop,
           "kappaBottom" = kappaBottom, "kappaTop" = kappaTop)
  }

  class(out) <- "clusterTree"
  out1 <- plotRule(out) ## relabel branches according to plotting rules.

  return(out1)
}
