# 02/07/11 -- Helper functions.

GoodnessOfFit <- function(x,type){
  ns <- NullModel(x,type=type)$n
  return(sum(((x-ns)^2)/ns, na.rm=TRUE))
}



ColorDendrogram <- function(hc, y, main = "", branchlength = 0.7, labels = NULL, xlab = NULL, sub = NULL, ylab = "", cex.main = NULL){
  if (is.null(labels))
  labels <- rep("", length(y))
  plot(hc, hang = 0.2, main = main, labels = labels, xlab = xlab,sub = sub, ylab = ylab, cex.main = cex.main)
  y <- y[hc$ord]
  if (is.numeric(y)) {
    y <- y + 1
    y[y == 7] <- "orange"
  }
  for (i in 1:length(hc$ord)) {
    o = hc$merge[, 1] == -hc$ord[i] | hc$merge[, 2] == -hc$ord[i]
    segments(i, hc$he[o] - branchlength, i, hc$he[o], col = y[i])
  }
}
        


permute.rows <- function(x){
  dd <- dim(x)
  n <- dd[1]
  p <- dd[2]
  mm <- runif(length(x)) + rep(seq(n) * 10, rep(p, n))
  matrix(t(x)[order(mm)], n, p, byrow = T)
}


balanced.folds <- function(y, nfolds = min(min(table(y)), 10)){
  totals <- table(y)
  fmax <- max(totals)
  nfolds <- min(nfolds, fmax)
  # makes no sense to have more folds than the max class size
  folds <- as.list(seq(nfolds))
  yids <- split(seq(y), y)
  # nice way to get the ids in a list, split by class
  ###Make a big matrix, with enough rows to get in all the folds per class
  bigmat <- matrix(NA, ceiling(fmax/nfolds) * nfolds, length(totals))
  for(i in seq(totals)) {
    bigmat[seq(totals[i]), i] <- sample(yids[[i]])
  }
  smallmat <- matrix(bigmat, nrow = nfolds) # reshape the matrix
  ### Now do a clever sort to mix up the NAs
  smallmat <- permute.rows(t(smallmat)) ### Now a clever unlisting
  x <- apply(smallmat, 2, function(x) x[!is.na(x)])
   if(is.matrix(x)){
    xlist <- list()
    for(i in 1:ncol(x)){
      xlist[[i]] <- x[,i]
    }
    return(xlist)
  }
  return(x)
}



se <- function(vec){
  return(sd(vec)/sqrt(length(vec)))
}


Soft <- function(x,a){
  return(sign(x)*pmax(abs(x)-a,0))
}

GetD <- function(ns, x, y, rho,beta,rhos=NULL){
  if(!is.null(rho) && !is.null(rhos)) stop("do you want to use rho or rhos in GetD function???")
  if(is.null(rhos)){
    uniq <- sort(unique(y))
    ds <- matrix(1, nrow=length(uniq), ncol=ncol(x))
    for(k in 1:length(uniq)){
      a <- colSums(x[y==uniq[k],])+beta
      b <- colSums(ns[y==uniq[k],])+beta
      ds[k,] <- 1+Soft(a/b-1,rho/sqrt(b))
    }
    return(ds)
  } else {
    uniq <- sort(unique(y))
    ds.list <- list()
    for(rho in rhos){
      ds <- matrix(1, nrow=length(uniq), ncol=ncol(x))
      for(k in 1:length(uniq)){
        a <- colSums(x[y==uniq[k],])+beta
        b <- colSums(ns[y==uniq[k],])+beta
        ds[k,] <- 1+Soft(a/b-1,rho/sqrt(b))
      }
      ds.list[[which(rhos==rho)]] <- ds
    }
    return(ds.list)
  }
}


print.poicla <- function(x,...){
  if(is.null(x$rhos) && is.null(x$rho)){
    if(!is.null(x[[1]]$alpha)) cat("Value of alpha used to transform data: ", x[[1]]$alpha, fill=TRUE)
    if(is.null(x[[1]]$alpha)) cat("No transformation performed.", fill=TRUE)
    cat("Type of normalization performed: ", x[[1]]$type, fill=TRUE)
    cat("Number of training observations: ", nrow(x[[1]]$x), fill=TRUE)
    if(!is.null(x[[1]]$xte)) cat("Number of test observations: ", nrow(x[[1]]$xte), fill=TRUE)
    cat("Number of features: ", ncol(x[[1]]$x), fill=TRUE)
    cat("-------------------------", fill=TRUE)
    cat("-------------------------", fill=TRUE)
    for(i in 1:length(x)){
      cat("-------------------------", fill=TRUE)
      cat("Value of rho used: ", x[[i]]$rho, fill=TRUE)
      cat("Number of features used in classifier: ", sum(apply(x[[i]]$ds!=1, 2, sum)!=0), fill=TRUE)
      if(sum(apply(x[[i]]$ds!=1, 2, sum)!=0)<100) cat("Indices of features used in classifier: ", which(apply(x[[i]]$ds!=1, 2, sum)!=0), fill=TRUE)
    }
  } else {
    if(!is.null(x$alpha)) cat("Value of alpha used to transform data: ", x$alpha, fill=TRUE)
    if(is.null(x$alpha)) cat("No transformation performed.", fill=TRUE)    
    cat("Type of normalization performed: ", x$type, fill=TRUE)    
    cat("Number of training observations: ", nrow(x$x), fill=TRUE)
    if(!is.null(x$xte)) cat("Number of test observations: ", nrow(x$xte), fill=TRUE)
    cat("Number of features: ", ncol(x$x), fill=TRUE)
    cat("Value of rho used: ", x$rho, fill=TRUE)
    cat("Number of features used in classifier: ", sum(apply(x$ds!=1, 2, sum)!=0), fill=TRUE)
    if(sum(apply(x$ds!=1, 2, sum)!=0)<100) cat("Indices of features used in classifier: ", which(apply(x$ds!=1, 2, sum)!=0), fill=TRUE)
  }
}
