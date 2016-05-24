
TSDatabaseDistances <- function(X, Y=NULL, distance, ...) {
  
#  If Y does not appear, redefine input parameters
 if (is.character(Y)) {
    distance <- Y
    Y <- NULL
  }

possible.distances <- c("euclidean", "manhattan", "minkowski", "infnorm",   
                        "ccor","sts", "dtw", "lb.keogh", "edr", "erp", "lcss", 
                        "fourier", "tquest", "dissim", "acf", "pacf", "ar.lpc.ceps",
                        "ar.mah", "ar.mah.statistic", "ar.mah.pvalue", "ar.pic",
                        "cdm", "cid", "cor", "cort", "wav", "int.per", "per", 
                        "mindist.sax", "ncd", "pred", "spec.glk", "spec.isd", 
                        "spec.llr", "pdc", "frechet")
  
  
distance <- match.arg(distance, possible.distances)
  
  # Initial checks
  if (! is.numeric(X) & ! is.matrix(X) & ! is.mts(X) &  ! is.zoo(X) & !
        is.xts(X) & ! is.list(X)) {
    stop('X must be a matrix, mts, zoo, xts or list object.')
  }  

  if (is.mts(X)) {
    X <- t(X)
    tx <- as.numeric(time(X))
  }
  if (is.zoo(X) | is.xts(X)) {
    X <- t(X)
    tx <- as.numeric(index(X))
  } else {
    tx<-NULL
  }
  
  if (! is.list(X)) {
    if (dim(X)[1] <= 1) {
      stop('The database must contain more than one series.')}
} else {
    if (length(X)<= 1) {
      stop('The database must contain more than one series.')}  
  }
  

  # Distance calculations for only one database
  if (is.null(Y)) {
    # Calculate distance matrix
    # Special cases of TSclust (more efficient than using dist.)
    if (distance ==  "ar.mah") {
      d1 <- dist(X, method="TSDistances", distance="ar.mah.statistic")
      d2 <- dist(X, method="TSDistances", distance="ar.mah.pvalue")
      d <- list(statistic=d1, pvalue=d2)
    } else if (distance ==  "ar.pic") {
      d <- as.dist(PairwiseDistances1(X, PairwiseARPicDistance, ...))
    } else if (distance ==  "ar.lpc.ceps") {
      d <- as.dist(PairwiseDistances1(X, PairwiseARLPCCepsDistance, ...))
    } else if (distance ==  "pred") {
      d <- as.dist(PairwisePredDistance(X, Y=NULL, ...))
    } else if (distance == "spec.llr") {
      d <- as.dist(PairwiseSpecLLRDistance(X, Y=NULL, ...))
    } else if (distance == "spec.glk") {
      d <- as.dist(PairwiseSpecGLKDistance(X, Y=NULL, ...))
    } else if (distance == "spec.isd") { 
      d <- as.dist(PairwiseSpecISDDistance(X, Y=NULL, ...))
    } else if (distance == "cdm") { 
      d <- as.dist(PairwiseDistances1(X, PairwiseCDMDistance, ...))
    } else if (distance == "ncd") { 
      d <- as.dist(PairwiseDistances1(X, PairwiseNCDDistance, ...))
    } else if (distance == "wav") {
      options(show.error.messages=FALSE)
      d <- diss.DWT(X)
      
    # For the PDC distance we use the original function from package pdc.
    # Faster than  using dist.
    } else if (distance == "pdc") {
      d <- pdcDist(t(X))
    } else if (distance == "frechet") {
      ty <- tx
      d <- as.dist(PairwiseDistances1(X, PairwiseFrechetDistance, ...))
      # For the rest of the cases: we use dist.
  } else {
      #options(show.error.messages=FALSE)
      d <- dist(X, method="TSDistances", distance=distance, tx=tx, ty=tx, ...)
    }
   
  # Calculate pairwise distances between series of 2 databases.
  # For TRAIN/TEST environments
} else {
    # Initial checks.
    if (! is.numeric(Y) & ! is.matrix(Y) & ! is.mts(Y) &  ! is.zoo(Y) & ! is.xts(Y) & ! is.list(Y)) {
      stop('Y must be a matrix, mts, zoo, xts or list object.')
    } 
    if (is.mts(Y)) {
      Y <- t(Y)
      ty <- as.numeric(time(Y))
    } else if (is.zoo(Y) | is.xts(Y)) {
      Y <- t(Y)
      ty <- as.numeric(index(Y))
  } else {
      ty <- NULL
    }
    
    if (! is.list(Y)) {
      if (dim(Y)[1] <= 1) {
        stop('The database must contain more than one series.')}
  } else {
      if (length(Y) <= 1) {
        stop('The database must contain more than one series.')}  
    }
    
    # Special cases of TSclust (more efficient than using dist.)
    if (distance == "ar.mah") {
      d1 <- dist(X, Y, method="TSDistances", distance="ar.mah.statistic")
      d2 <- dist(X, Y, method="TSDistances", distance="ar.mah.pvalue")
      d <- list(statistic=d1, pvalue=d2)
    } else if (distance == "ar.pic") {
      d <- PairwiseDistances2(X, Y, PairwiseARPicDistance, ...)
    } else if (distance == "ar.lpc.ceps") {
      d <- PairwiseDistances2(X, Y, PairwiseARLPCCepsDistance, ...)
    } else if (distance == "pred") {
      d <- PairwisePredDistance(X, Y, ...)
    } else if (distance == "spec.llr") {
      d <- PairwiseSpecLLRDistance(X, Y, ...)
    } else if (distance == "spec.glk") {
      d <- PairwiseSpecGLKDistance(X, Y, ...)
    } else if (distance == "spec.isd") { 
      d <- PairwiseSpecISDDistance(X, Y, ...)
    } else if (distance == "cdm") { 
      d <- PairwiseDistances2(X, Y, PairwiseCDMDistance, ...)
    } else if (distance == "ncd") { 
      d <- PairwiseDistances2(X, Y, PairwiseNCDDistance, ...)
    } else if (distance == "frechet") {
      d <- PairwiseDistances2(X, Y, PairwiseFrechetDistance, ...)
      # For the rest of the cases: we use dist.
    } else if (distance == "pdc") {
      d <- PDCDist2(t(X), t(Y))
    # Both the training and testing databases are used for feature extraction.
    } else if (distance == "wav") {
      if (is.list(X) & is.list(Y)) {
        X <- matrix(unlist(X), ncol = length(X[[1]]), byrow = TRUE)
        Y <- matrix(unlist(Y), ncol = length(Y[[1]]), byrow = TRUE)  
    }
      d <- as.matrix(diss.DWT(rbind(X, Y)))
    
      n <- dim(X)[1]
      m <- dim(Y)[1]
      d <- d[-c(1:n),-c(1:m)]
      # Rest of the cases: we use dist.
  } else {
      #options(show.error.messages=FALSE)
      d <- dist(X, Y, method="TSDistances", distance=distance, tx=tx, ty=ty, ...)
    }
  }
  options(show.error.messages=TRUE)
  return(d) 
}
 
   
  