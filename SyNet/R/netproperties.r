netproperties <- function(mt, cutoff = 0, dichotomization = ">") {
  stopifnot(is.matrix(mt))
  stopifnot(colnames(mt) == rownames(mt))
  stopifnot(isSymmetric(mt))
  nsp <- nrow(mt)
  mt <- do.call(dichotomization, list(mt, cutoff))
  stopifnot(is.logical(mt))
  diag(mt) <- TRUE
  adj <- apply(mt, 1, which)
  if(!is.list(adj)) adj <- data.frame(matrix(adj, ncol = nsp))
  #Identify components
  NodeComponent <- rep(0, nsp)
  maxc <- 1
  dfs <- function (i) {
      NodeComponent[i] <<- maxc
      for(j in adj[[i]]) if(NodeComponent[j] == 0) dfs(j)
  }
  while(1) {
      root <- match(0, NodeComponent, nomatch = 0)
      if(root==0) break
      dfs(root)
      maxc <- max(NodeComponent) + 1
  }
  #Identify cliques
  maxcliq <- function() {
      # The total number of vertices is provided by nsp variable
      MC <- c() # storage for maximal cliques
      R <- c() # currently growing clique
      P <- 1:nsp # prospective nodes connected to all nodes in R
      X <- c() # nodes already processed
      BKv2 <- function (R, P, X){
          if (length(P)==0 && length(X) == 0) {
              # report R as a maximal clique
              newMC <- rep(0, nsp)
              newMC[R] <- 1 # newMC contains ones at indices equal to the values in R
              MC <<- cbind(MC, newMC)}
          else {
              # choose pivot
              ppivots = union(P, X) # potential pivots
              binP = rep(0, nsp)
              binP[P] <- 1 # binP contains ones at indices equal to the values in P
              # rows of mt(ppivots,:) contain ones at the neighbors of ppivots
              pcounts <- array(mt[ppivots,,drop = FALSE]%*%binP)  
              # cardinalities of the sets of neighbors of each ppivots intersected with P
              ind = which.max(pcounts)
              u_p = ppivots[ind] # select one of the ppivots with the largest count
              for (u in intersect(which(!mt[u_p,]), P)) { 
                  # all prospective nodes who are not neighbors of the pivot
                  P <- setdiff(P, u)
                  Rnew <- c(R, u)
                  Nu <- which(mt[u,])
                  Pnew <- intersect(P, Nu)
                  Xnew <- intersect(X, Nu)
                  BKv2(Rnew, Pnew, Xnew)
                  X <- c(X, u)
              }
          }
      }
      BKv2(R, P, X)
      return(MC)
  }
  diag(mt) <- FALSE
  Cliques <- maxcliq()
  colnames(Cliques) <- paste("Clique_", 1:ncol(Cliques), sep = "")
  #Identify geodesic distances between nodes
  APD <- function(mt, n = nrow(mt)) {
     Z <- mt%*%mt
     degree <- diag(Z)
     B <- ifelse((mt + Z) > 0, 1, 0)
     diag(B) <- 0
     if (all(B[lower.tri(B)]==1)) return (2*B-mt)
     T <- APD(B)
     X <- T%*%mt
     D <- 2*T
     for (i in 1:(n -1))
       for (j in (i + 1):n)
       if (X[i,j] < T[i,j]*degree[j]) D[i,j] <- D[j,i] <- D[i,j] - 1
     diag(D) <- 0
     return (D)
  }
  diag(mt) <- 0
  AllPairsDistance <- matrix(NA, nsp, nsp)
  diag(AllPairsDistance) <- 0
  for(comp in 1:max(NodeComponent)) {
    whsp <- which(NodeComponent == comp)
    if(length(whsp) == 1) next
    AllPairsDistance[whsp, whsp] <- APD(mt[whsp, whsp])
  }
  #Calculate the degrees of nodes
  NodeDegree <- unlist(lapply(adj, length)) - 1
  out <- c()
  out$Adjacency <- mt
  out$Components <- NodeComponent
  out$Degree <- NodeDegree
  out$Geodesic <- AllPairsDistance
  out$Cliques <- Cliques     
  return(out)
}



