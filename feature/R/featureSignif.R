
featureSignif <- function(x, bw, gridsize, scaleData=FALSE, addSignifGrad=TRUE, addSignifCurv=TRUE, signifLevel=0.05)                  
{ 
  ## tau is effective kernel support
  tau <- 5
  if (is.vector(x))
  {
    d <- 1
    n <- length(x)
    names.x <- deparse(substitute(x))
    if (scaleData)  x <- (x-min(x))/(max(x) - min(x))
  }
  else
  {  
    d <- ncol(x)
    n <- nrow(x)
    names.x <- colnames(x)
    
    if (is.null(names.x))
    {
      names.xx <- deparse(substitute(x))
      names.xx <- strsplit(names.xx, "\\[")[[1]][1]
      names.x <- paste(names.xx, "[,", 1:d,"]",sep="") 
    }
    if (scaleData)
      for (i in 1:d)
        x[,i] <- (x[,i]-min(x[,i]))/(max(x[,i]) - min(x[,i]))
  }
  x <- as.matrix(x)
  if (d>4)
    stop("Feature significance currently only available for 1- to 4-dimensional data")
  
  if (missing(gridsize)) 
  {
    if (d==1) gridsize <- 401
    if (d==2) gridsize <- rep(151,2)
    if (d==3) gridsize <- rep(51,3)
    if (d==4) gridsize <- rep(21,4)
  }
  
  ## Set some defaults

  if (missing(bw))           ## b/w not specified -> interactive 
  {
    bw.range <- dfltBWrange(as.matrix(x),tau)
    bw <- matrix(unlist(bw.range), nrow=2, byrow=FALSE)
    dfltCounts.out <- dfltCounts(x, gridsize, apply(bw, 2, max))
   
    h.low <- bw[1,]
    h.upp <- bw[2,]
    hmix.prop <- 1/4
    h.init <- h.low^(hmix.prop)*h.upp^(1-hmix.prop)
    h <- h.init
  }
  else
  {
    dfltCounts.out <- dfltCounts(x,gridsize, bw)
    h <- bw
  }
  gcounts <- dfltCounts.out$counts
  ##range.x <- lapply(dfltCounts.out$eval.points, range)
  range.x <- dfltCounts.out$range.x

  dest <- drvkde(gcounts, drv=rep(0,d), bandwidth=h, binned=TRUE, range.x=range.x, se=FALSE, gridsize=gridsize)
  dest$est[dest$est<0] <- 0 
 
  ## significant features 
  
  SignifFeatureRegion.mat <- SignifFeatureRegion(n,d,gcounts,gridsize,dest, h, signifLevel, range.x, grad=addSignifGrad, curv=addSignifCurv)
  ESS <- n*dest$est*prod(h)*(sqrt(2*pi)^d)
  SigESS <- ESS >= 5

  SignifGradRegion.mat <- SignifFeatureRegion.mat$grad & SigESS
  SignifGradData.mat <- SignifFeatureData(x, d, dest,SignifGradRegion.mat)
  SignifGradDataPoints <- x[SignifGradData.mat,]
  
  SignifCurvRegion.mat <- SignifFeatureRegion.mat$curv & SigESS
  SignifCurvData.mat <- SignifFeatureData(x, d,  dest,SignifCurvRegion.mat) 
  SignifCurvDataPoints <- x[SignifCurvData.mat,] 
    
  feat <- c(list(x=x, names=names.x, bw=h, fhat=dest), SignifFeatureRegion.mat, list(gradData=SignifGradData.mat, gradDataPoints=SignifGradDataPoints, curvData=SignifCurvData.mat, curvDataPoints=SignifCurvDataPoints))
    
  class(feat) <- "fs"

  return(feat)
}
  



