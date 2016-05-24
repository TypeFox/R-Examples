piecewiselin <- function(distmatrix, maxdist=0.1*max(distmatrix)){
  ncd <- ncol(distmatrix)
  nmatrix <- distmatrix/maxdist  
  for (i in 1:(ncd-1))
    for (j in (i+1):ncd)
      if (nmatrix[i,j]>1)
        nmatrix[i,j] <- nmatrix[j,i] <- 1
  nmatrix
}
          
# Geographic Kulczynski distance, rows are regions
geco <- function(regmat,
                   geodist=as.dist(matrix(as.integer(!diag(nrow(regmat))))),
                   transform="piece",
                   tf=0.1,
                   countmode=ncol(regmat)+1){
#  print(tf)
  nart <- ncol(regmat)
  ncell <- nrow(regmat)
  jdist <- rep(0, nart * nart)
  dim(jdist) <- c(nart, nart)
  if (transform=="none")
    geodist <- as.matrix(geodist)
  if (transform=="log")
    geodist <- log(as.matrix(tf*geodist)+1)
  if (transform=="sqrt")
    geodist <- sqrt(as.matrix(tf*geodist))
  if (transform=="piece"){
    mg <- max(geodist)
    geodist <- piecewiselin(as.matrix(geodist),tf*mg)
  }
  fi <- list()
  nci <- numeric(0)
  for (i in 1:nart){
    fi[[i]] <- (1:ncell)[as.logical(regmat[,i])]
    nci[[i]] <- sum(regmat[, i])
  }
  for (i in 1:(nart - 1)) {
    if (round(i/countmode)==(i/countmode))
      cat("Computing gb-distances for  species ",i,"\n")
    for (j in (i + 1):nart) {
      nsi <- numeric(0)
      gfi <- geodist[fi[[i]],fi[[j]],drop=FALSE]
#      print(gfi)
      nsi[1] <- sum(apply(gfi,1,min))
      nsi[2] <- sum(apply(gfi,2,min))
      jdist[i,j] <- jdist[j,i] <- (nsi[1]/nci[i]+nsi[2]/nci[j])/2
      if (is.na(jdist[i, j])) 
        cat("Warning! NA at i=", i, ", j=", j, "\n")
    }
  }
  jdist
}

geo2neighbor <- function(geodist,cut=0.1*max(geodist)){
  geodist <- as.matrix(geodist)
  n <- nrow(geodist)
  nblist <- list()
  for (i in 1:n) nblist[[i]] <- numeric(0)
  for (i in 1:(n-1))
    for(j in (i+1):n)
      if (geodist[i,j]<=cut){
        nblist[[i]] <- c(nblist[[i]],j)
        nblist[[j]] <- c(nblist[[j]],i)
      }
  nblist
}
        

hprabclust <- function (prabobj, cutdist=0.4, cutout=1,
                        method="average", nnout=2, mdsplot=TRUE,
                        mdsmethod="classical") 
{
    cf <- match.call()
#    if (mdsplot & mdsmethod!="classical")
#      require(MASS)
    # "Data-alphabetical" ordering 
    oregions <- order(prabobj$specperreg)
    prabo1 <- prabobj$prab[oregions,]
    ospecies <- do.call("order",as.data.frame(t(prabo1)))    
    dma <- prabobj$distmat[ospecies,ospecies]
    if (mdsplot){
      mdsdim=2
      if (mdsmethod != "classical") {
          mindm <- min(dma[dma > 0])/10
          for (i in 1:(prabobj$n.species - 1))
            for (j in (i + 1):prabobj$n.species) if (dma[i, j] < mindm) 
              dma[i, j] <- dma[j, i] <- mindm
      }
      mdsout <-
        switch(mdsmethod, classical = cmdscale(dma, k = mdsdim), 
        kruskal = isoMDS(dma, k = mdsdim), sammon = sammon(dma, 
            k = mdsdim))
      if (mdsmethod=="classical") mds <- mdsout
      else mds <- mdsout$points
    }
    else mds <- NULL
    n <- prabobj$n.species
    nnd <- c()
    nout <- rep(TRUE,n)
    for (i in 1:n){
      nnd[i] <- sort(dma[i, ])[nnout + 1]
      if (nnd[i]>cutout) nout[i] <- FALSE
    }
    noisen <- n-sum(nout)
    dm <- as.dist(dma[nout,nout])
    cl1 <- hclust(dm, method=method)
    rclustering <- cl2 <- cutree(cl1, h=cutdist)
    nc <- max(cl2)
#    nout[ospecies] <- nout
#    ncl <- max(cl2)
    csum <- function(nx, cv) {
        out <- c()
        for (i in 1:length(nx)) out[i] <- sum(cv == nx[i])
        out
    }
    cs <- csum(1:nc, cl2)
    ocs <- order(-cs)
    for (i in 1:nc) cl2[rclustering == ocs[i]] <- i
    clustering <- rep(0,n)
    clustering[ospecies][nout] <- cl2
    nmr <- max((1:nc)[cs[ocs]>nnout])
#    print(cs[ocs])
#    print(nmr)
    rclustering <- rep(0,n)
    rclustering[clustering<=nmr] <- clustering[clustering<=nmr]
    symbols <- c("N")
    if (nmr>0) symbols <- c("N", sapply(1:nmr, toString))
    clsym <- symbols[rclustering+1]
    if (mdsplot){
      mds[ospecies,] <- mds
      plot(mds[,1:2],pch=clsym)
    }
    out <- list(clustering = clustering, rclustering=rclustering,
                cutdist=cutdist, method=method,
                cutout=cutout,nnout=nnout,noisen=noisen,
                symbols = clsym, points=mds, call=cf)
    class(out) <- "comprabclust"
    out
}

"print.comprabclust" <-
function(x, ...){
  cat("* Clustered presence-absence matrix * \n\n")
  cat("Clustered by hclust with noise detection, method=: ", x$method, "\n\n")
  cat("Distance value to cut tree:", x$cutdist,"\n")
  cat("Minimum cluster size (below is noise):", x$nnout+1,"\n")
  cat("Call: \n")
  print(x$call)
  cat("\n Cluster memberships:\n")
  print(x$rclustering)
  cat("\n")
}
