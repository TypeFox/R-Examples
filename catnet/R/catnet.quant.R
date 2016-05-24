########################################################################
# Categorical Network Class Methods
# Discretization

cnDiscretize <- function(data, numCategories, mode = "uniform", qlevels = NULL) {

  if(!is.matrix(data) && !is.data.frame(data))
    stop("data should be a matrix or data frame")

  if(is.matrix(data)) {
    samples <- data
    outdataframe <- FALSE
  }
  else {
    samples <- t(data)
    outdataframe <- TRUE
  }
  
  numnodes <- dim(samples)[1]
  numsamples <- dim(samples)[2]
  if(numnodes < 1 || numsamples < 1)
    stop("Invalid sample")
  
  if(is.vector(numCategories)) {
    len <- length(numCategories)
    numcats <- floor(numCategories[len])
    numcats <- rep(numcats, numnodes)
    if(len <= numnodes) {
      for(j in 1:len)
        numcats[j] <- numCategories[j]
    }
    else {
      for(j in 1:numnodes)
        numcats[j] <- numCategories[j]
    }
    for(j in 1:numnodes)
        if(numcats[j] < 2)
          numcats[j] <- 2
  }
  else {
    numcats <- floor(numCategories)
    if(numcats < 2) {
      warning("set numCategories to 2")
      numcats <- 2
    }
    numcats <- rep(numcats, numnodes)
  }

  if(mode == "quantile") {
    if(is.null(qlevels) || !is.list(qlevels) || length(qlevels) < numnodes) {
      qlevels <- vector("list", numnodes)
      for(j in 1:numnodes)
        qlevels[[j]] <- seq(1:(numcats[j]-1))/numcats[j]
    }
    for(j in 1:numnodes) {
      if(length(qlevels[[j]]) < numcats[j]-1 || min(qlevels[[j]]) < 0 || max(qlevels[[j]]) > 1)
        qlevels[[j]] <- seq(1:(numcats[j]-1))/numcats[j]
    }
  }
  if(mode == "uniform") {
    if(is.null(qlevels) || !is.list(qlevels) || length(qlevels) < numnodes) {
      qlevels <- vector("list", numnodes)
      for(j in 1:numnodes)
        qlevels[[j]] <- rep(1, numcats[j])
    }
    for(j in 1:numnodes) {
      sl <- sum(qlevels[[j]])
      if(length(qlevels[[j]]) < numcats[j] || sl <= 0) {
        qlevels[[j]] <- rep(1, numcats[j])
        sl <- sum(qlevels[[j]])
      }
      qlevels[[j]] <- qlevels[[j]]/sl
    }
  }

  data <- matrix(rep(NA, length(samples)), nrow=dim(samples)[1])
  for(j in 1:numnodes) {
    col <- samples[j,]
    ids <- !is.na(col)
    col <- as.numeric(col)
    if(!is.numeric(col))
      stop("data should be numeric")
    qq <- qlevels[[j]]
    if(mode == "quantile") {
      qq <- quantile(col[ids], qlevels[[j]])
    }
    if(mode == "uniform") {
      if(sum(ids) > 0) {
        minc <- min(col[ids])
        maxc <- max(col[ids])
        qq <- sapply(1:numcats[j], function(i) minc+(maxc-minc)*sum(qlevels[[j]][1:i]))
      }
    }
    for(i in 1:length(col)) {
      if(!ids[i]) 
        next
      id <- which(qq > col[i])
      if(length(id)>0)
        data[j,i] <- min(id)
      else
        data[j,i] <- numcats[j]
    }
    data[j,] <- as.integer(data[j,])
  }

  data <- matrix(as.integer(data), nrow=numnodes)
  rownames(data) <- rownames(samples)
  colnames(data) <- colnames(samples)

  if(outdataframe) {
    data <- as.data.frame(t(data))
    ## R converts integer back to numeric, hate it
    for(j in 1:numnodes)
      data[,j] <- data[,j]
  }

  return(data)
}

cnGMM <- function(data, numCategories, maxiter=5, equal.sigma=FALSE) {

  if(!is.vector(data) && !is.numeric(data))
    stop("data should be a numeric vector")

  n <- length(data)

  if(numCategories<1)
     numCategories <- 1
  cc <- numCategories
  if(maxiter < 2)
    maxiter <- 2
  
  qq <- c(min(data), quantile(data, (1:(cc-1))/cc), max(data))
  names(qq) <- 1:(cc+1)
  mu <- (qq[1:cc]+qq[2:(cc+1)])/2
  sig <- ((qq[1:cc]-qq[2:(cc+1)])/2)^2
  w <- rep(1,cc)/cc
  wc <- matrix(rep(0,cc*n),nrow=cc)

  if(maxiter < 1) {
    if(equal.sigma)
      sig <- rep(mean(sig), cc)
    return(list(w=w, mu=mu, sigma=sqrt(sig)))
  }

  eps <- 1e-8
  for(iter in 1:maxiter) {
    oldmu <- mu
    ##cat(mu,"\n")
    ##cat(iter,": ", sig,"\n")
    if(equal.sigma)
      sig <- rep(mean(sig), cc)
    for(c in 1:cc) 
      wc[c,] <- w[c]*exp(-(data-mu[c])^2/sig[c])/sqrt(sig[c])
    for(j in 1:n) 
      wc[,j] <- wc[,j]/sum(wc[,j])
    for(c in 1:cc) {
      mu[c] <- sum(wc[c,]*data)/sum(wc[c,])
      ##sig[c] <- sum(wc[c,]*data^2)/sum(wc[c,])-mu[c]^2
      sig[c] <- sum(wc[c,]*(data-mu[c])^2)/sum(wc[c,])
      if(sig[c] < eps)
        sig[c] <- eps
      w[c] <- mean(wc[c,])
    }
    if(sum((oldmu-mu)^2) < eps)
      break
  }
  if(equal.sigma)
      sig <- rep(mean(sig), cc)
  return(list(w=w, mu=mu, sigma=sqrt(sig)))
}

cnUMM <- function(data, numCategories) {

  if(!is.vector(data) && !is.numeric(data))
    stop("data should be a numeric vector")

  n <- length(data)
    
  if(numCategories<1)
     numCategories <- 1
  cc <- numCategories
  if(maxiter < 2)
    maxiter <- 2
  
  qq <- min(data) + (max(data)-min(data))*(0:cc)/cc
  names(qq) <- 1:(cc+1)
  mu <- (qq[1:cc]+qq[2:(cc+1)])/2
  sig <- ((qq[1:cc]-qq[2:(cc+1)])/2)^2
  w <- rep(1,cc)/cc

  qq[cc+1] <- qq[cc+1]+1
  for(i in 1:cc) {
    ind <- data>=qq[i] & data<qq[i+1]
    w[i] <- sum(ind)
  }
  w <- w/sum(w)
  
  return(list(w=w, mu=mu, sigma=sqrt(sig)))
}

cnQMM <- function(data, numCategories) {

  if(!is.vector(data) && !is.numeric(data))
    stop("data should be a numeric vector")

  n <- length(data)
    
  if(numCategories<1)
     numCategories <- 1
  cc <- numCategories
  if(maxiter < 2)
    maxiter <- 2

  qq <- c(min(data), quantile(data, (1:(cc-1))/cc), max(data))
  names(qq) <- 1:(cc+1)
  mu <- (qq[1:cc]+qq[2:(cc+1)])/2
  sig <- ((qq[1:cc]-qq[2:(cc+1)])/2)^2
  w <- rep(1,cc)/cc

  qq[cc+1] <- qq[cc+1]+1
  for(i in 1:cc) {
    ind <- data>=qq[i] & data<qq[i+1]
    w[i] <- sum(ind)
  }
  w <- w/sum(w)
  
  return(list(w=w, mu=mu, sigma=sqrt(sig)))
}

