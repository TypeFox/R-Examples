medianHL <- function(x, na.rm=FALSE, ...) UseMethod("medianHL")

medianHL.default <- function(x, na.rm=FALSE, method=c("HL1","HL2","HL3"), prop=NULL) {
  method <- match.arg(method)
  if (!is.null(prop))
    if (prop <= 0 | prop >=1)
      stop("'prop' is outside (0,1)")
  if (na.rm)
    x <- x[!is.na(x)]
  if ((n <- length(x))==0) {
    warning("No observations (at least after removing missing values)")
    return(NULL)
  }
  if (method=="HL2") {
    nt <- n*(n+1)/2
    nm <- n
  } else if (method=="HL1") {
    nt <- n*(n-1)/2
    nm <- n-1
  } else {
    nt <- n^2
    nm <- n
    im <- 1
  }
  if (is.null(prop)) {
    meanpairs <- rep(0,nt)
    ni <- 0
    for (i in 1:nm) {
      if (method=="HL1")
        im <- i+1
      else if (method=="HL2")
        im <- i
      for (j in im:n) {
        ni <- ni + 1
        meanpairs[ni] <- (x[i]+x[j])/2
      }
    }    
  } else {
    np <- round(nt*prop)
    indici <- sample(x=nt, size=np, replace=FALSE)
    if (np < 1)
      np <- 1
    meanpairs <- rep(0,np)
    ni <- 0
    npi <- 0
    for (i in 1:nm) {
      if (method=="HL1")
        im <- i+1
      else if (method=="HL2")
        im <- i
      for (j in im:n) {
        ni <- ni + 1
        if (any(indici==ni)) {
          npi <- npi + 1
          meanpairs[npi] <- (x[i]+x[j])/2
        }
      }
    }    
  }
  median.default(meanpairs)
}
