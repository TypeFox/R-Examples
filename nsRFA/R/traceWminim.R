

# ----------------------------------------------------------------------------- #

traceWminim <- function (X,centers) {

  # B.Everitt [1974] Cluster Analysis, pages 42-43

  # INPUT
  # X        A numeric matrix of data, or an object that can be coerced to
  #          such a matrix (such as a numeric vector or a data frame with
  #          all numeric columns).
  # centers  The number of clusters.

  X <- data.frame(X)    # to be sure that row.names exist  (12-10-2007)

  n <- dim(X)[1]
  k <- dim(X)[2]

  #X.norm <- (X - matrix(mean(X),nrow=n,ncol=k,byrow=TRUE))/matrix(sd(X),nrow=n,ncol=k,byrow=TRUE)

  d <- dist(X, method = "euclidean")
  tree <- hclust(d, method = "ward")

  clusters <- cutree(tree,centers)
  names(clusters) <- row.names(X)    # to be sure that row.names exist  (12-10-2007)

  fine=FALSE
  cont=1
  while(fine == FALSE) {
    traceW.1 <- sumtraceW(clusters,X)
    scambi <- nearest(clusters,X)
    traceW.2 <- rep(NA,centers)
    for (i in 1:centers) {
      clusters.mod <- clusters; clusters.mod[scambi[i]] <- i
      traceW.2[i] <- sumtraceW(clusters.mod,X)
    }
    min.traceW.2 <- min(traceW.2); pos.min.traceW.2 <- which(traceW.2==min(traceW.2))
    if(min.traceW.2 > traceW.1) {
      fine=TRUE
    }
    else {
      clusters[scambi[pos.min.traceW.2]] <- pos.min.traceW.2
    }
    cont <- cont+1
  }


  return(clusters)

}


# ------------------------------------------------------------------------------- #

sumtraceW <- function (clusters,X) {

  # INPUT
  # X         A numeric matrix of data, or an object that can be coerced to
  #           such a matrix (such as a numeric vector or a data frame with
  #           all numeric columns).
  # clusters  A numeric vector containing the subdivision of X in clusters

  #clusters <- round(as.numeric(clusters))
  k <- max(clusters)
  traceW <- rep(NA,k)
  for (i in 1:k) {
    gruppo <- X[clusters==i,]
    mean.vector <- apply(gruppo,2,mean)
    distanze <- dist(rbind(mean.vector,gruppo), method = "euclidean")
    distanze <- as.matrix(distanze)[,1][-1]
    W <- crossprod(t(distanze))
    traceW[i] <- sum(diag(W))
  }

  sumtraceW <- sum(traceW)

  return(sumtraceW)
}


# ------------------------------------------------------------------------------- #

nearest <- function (clusters,X) {

  # INPUT
  # X         A numeric matrix of data, or an object that can be coerced to
  #           such a matrix (such as a numeric vector or a data frame with
  #           all numeric columns).
  # clusters  A numeric vector containing the subdivision of X in clusters

  k <- max(clusters)
  near <- rep(NA,k)
  for (i in 1:k) {
    gruppo <- X[which(clusters==i),]
    altri <- X[-which(clusters==i),]
    mean.vector <- apply(gruppo,2,mean)
    distanze <- dist(rbind(mean.vector,altri), method = "euclidean")
    distanze <- as.matrix(distanze)[,1][-1]
    near[i] <- names(which(distanze==min(distanze)))
  }

  return(near)
}

