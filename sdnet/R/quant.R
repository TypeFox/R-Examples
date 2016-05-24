########################################################################
# Categorical Network Class Methods
# Discretization

get.range <- function(col, cover) {
  col <- as.numeric(col)
  if(!is.numeric(col))
    stop("data should be numeric")
  ind <- !is.na(col)
  if(cover > 0 && cover < 1) {
    scol <- sort(col)
    i1 <- 1
    i2 <- length(scol)
    q1 <- scol[i1]
    q2 <- scol[i2]
    while(i1 < i2 && (i2-i1)>cover*length(col)) {
      mm <- median(scol[i1:i2])
      if(abs(q1-mm) < abs(q2-mm))
        i2 <- i2-1
      else
        i1 <- i1+1
      q1 <- scol[i1]
      q2 <- scol[i2]
    }
    ##q1 <- quantile(col[ind], (1-cover)/2)
    ##q2 <- quantile(col[ind], 1-(1-cover)/2)
    ind <- ind & col>=q1 & col<= q2
  }
  return(ind)
}

hardDiscretize <- function(data, numcats=3, marginal = "uniform", learnset=NULL, cover=0.95, qlevels=NULL) {

  if(!is.matrix(data) && !is.data.frame(data))
    stop("data should be a matrix or data frame")

  if(is.matrix(data)) {
    samples <- data
    outdataframe <- FALSE
  }
  if(is.data.frame(data)) {
    samples <- t(data)
    outdataframe <- TRUE
  }
  
  numnodes <- dim(samples)[1]
  numSamples <- dim(samples)[2]
  if(numnodes < 1 || numSamples < 1)
    stop("No valid sample is specified.")

  if(!is.null(learnset)) {
    if(!is.integer(learnset))
      stop("learnset should be sample indices")
    if(length(learnset)>numSamples)
      learnset <- learnset[1:numSamples]
    if(sum(learnset<0)>0 || sum(learnset>numSamples)>0)
      stop("learnset should be sample indices")
  }
  else
    learnset <- 1:numSamples
  ## translation and scaling
  #if(marginal == "uniform" && length(learnset) < numSamples) {
  #  for(j in 1:numnodes) {
  #    col <- samples[j,learnset]
  #    ran <- range(col[get.range(col, cover)])
  #    if(ran[2]>ran[1])
  #      samples[j,learnset] <- 2 * (samples[j,learnset] - (ran[2]+ran[1])/2)/ (ran[2]-ran[1])
  #    col <- samples[j,-learnset]
  #    ran <- range(col[get.range(col, cover)])
  #    if(ran[2]>ran[1])
  #      samples[j,-learnset] <- 2 * (samples[j,-learnset] - (ran[2]+ran[1])/2)/ (ran[2]-ran[1])
  #  }
  #}
 
  if(is.vector(numcats)) {
    len <- length(numcats)
    numcats <- floor(numcats[len])
    numcats <- rep(numcats, numnodes)
    if(len <= numnodes) {
      for(j in 1:len)
        numcats[j] <- numcats[j]
    }
    else {
      for(j in 1:numnodes)
        numcats[j] <- numcats[j]
    }
    for(j in 1:numnodes)
        if(numcats[j] < 2)
          numcats[j] <- 2
  }
  else {
    numcats <- floor(numcats)
    if(numcats < 2) {
      warning("set numcats to 2")
      numcats <- 2
    }
    numcats <- rep(numcats, numnodes)
  }

  if(marginal == "quantile") {
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
  if(marginal == "uniform") {
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
    if(!is.null(learnset))
      col <- samples[j,learnset]
    ids <- !is.na(col)
    if(sum(ids)<2) {
      stop("not enough data available: ", col)
    }
    col <- as.numeric(col)
    if(!is.numeric(col))
      stop("data should be numeric")
    ins <- !is.na(col)
    qq <- qlevels[[j]]
    if(marginal == "quantile") {
      qq <- quantile(col[ids], qlevels[[j]])
    }
    if(marginal == "uniform") {
      ids <- get.range(col, cover)
      if(sum(ids) > 0) {
        minc <- min(col[ids])
        maxc <- max(col[ids])
        qq <- sapply(1:numcats[j], function(i) minc+(maxc-minc)*sum(qlevels[[j]][1:i]))
      }
    }
    if(!is.null(learnset)) {
      col <- samples[j,]
      ids <- !is.na(col)
      col <- as.numeric(col)
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

softDiscretize <- function(data, numcats=3, marginal = "uniform", learnset=NULL, cover=0.95, maxiter=100, eps=1e-15, weights=NULL) {

  if(!is.matrix(data) && !is.data.frame(data))
    stop("data should be a matrix or data frame")
  if(is.matrix(data)) {
    outdataframe <- FALSE
  }
  if(is.data.frame(data)) {
    data <- t(data)
    outdataframe <- TRUE
  }
  if(!is.numeric(data))
    stop("data should be numeric")

  numnodes <- dim(data)[1]
  numSamples <- dim(data)[2]
  if(numnodes < 1 || numSamples < 1)
    stop("No valid sample is specified.")

  if(!is.null(learnset)) {
    if(!is.integer(learnset))
      stop("learnset should be sample indices")
    if(length(learnset)>numSamples)
      learnset <- learnset[1:numSamples]
    if(sum(learnset<0)>0 || sum(learnset>numSamples)>0)
      stop("learnset should be sample indices")
  }
  if(is.null(learnset))
    learnset <- 1:numSamples

  maxcats <- max(numcats)
  if(maxcats < 0)
    stop("Invalid numcats parameter")
  if(!is.vector(numcats) || (is.vector(numcats) && length(numcats) != numnodes))
    numcats <- rep(maxcats, numnodes)

  maxiter <- as.integer(maxiter)

  if(is.null(weights) || sum(weights) <= 0)
    weights <- rep(1/numSamples, numSamples)
  weights[weights < 0] <- 0
  weights <- weights / sum(weights)

  if(marginal == "gauss" || marginal == "uniform" || marginal == "quantile") {
    plist <- .Call("ccnSoftQuant", 
                          data, weights, numcats, learnset, cover, marginal, maxiter, eps, 
                          PACKAGE="sdnet")
    pdata <- plist[[1]]
    ddata <- plist[[2]]
    pdata <- matrix(pdata, nrow=numcats*numnodes)
    ddata <- matrix(ddata, nrow=numnodes)
  
    colnames(ddata) <- colnames(data)
    rownames(ddata) <- rownames(data)
    colnames(pdata) <- colnames(data)

    if(outdataframe) {
      ddata <- as.data.frame(t(ddata))
      pdata <- as.data.frame(t(pdata))
      ## R converts integer back to numeric, hate it
      for(j in 1:numnodes)
        ddata[,j] <- ddata[,j]
    }
    return(list(ddata=ddata, pdata=pdata, mu=t(matrix(plist[[3]],nrow=numcats)), sig=t(matrix(plist[[4]],nrow=numcats)), w=t(matrix(plist[[5]],nrow=numcats))))
  }
  return(NULL)
}

cnDiscretize <- function(data, numcats=3, mode="soft", marginal="quantile", learnset=NULL, cover=0.95, maxiter=100, eps=1e-8, weights=NULL) {
  if(marginal != "uniform" && marginal != "quantile" && marginal != "gauss")
    stop("Unsupported marginal")
  if(mode == "hard")
    return(hardDiscretize(data, numcats, marginal = marginal, learnset=learnset, cover=cover, qlevels=NULL))
  if(mode == "soft")
    return(softDiscretize(data, numcats, marginal = marginal, learnset=learnset, cover=cover, maxiter=maxiter, eps, weights))
  stop("Unsupported mode")
}

