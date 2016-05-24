runmbo <- function(x, y, fr = 1, est = "mom", nboot = 40){
                  
  #
  # running interval smooth with bagging
  #
  # fr controls amount of smoothing
  # tr is the amount of trimming
  #
  # Missing values are automatically removed.
  #
  # RNA=F, do not remove missing values when averaging
  # (computing the smooth) at x
  # xout=T removes points for which x is an outlier
  # eout=F removes points for which (x,y) is an outlier
  # nmin  estimate y|x only when number of points close
  # to x is > nmin
  # atr is amount of trimming when averaging over the bagged
  # values
  # est is the measure of location to be estimated
  # est=tmean means estimate 20% trimmed mean of y given x
  #
  est <- match.arg(est, c("mom", "onestep", "median"), several.ok = FALSE)
  est <- get(est)                 
  temp <- cbind(x,y)
  temp <- elimna(temp) # Eliminate any rows with missing values
  #if(eout && xout)stop("Not allowed to have eout=xout=T")
  #if(eout){
  #  flag<-outfun(temp,plotit=FALSE)$keep
  #  temp<-temp[flag,]
  #}
  #if(xout){
  #  flag<-outfun(x,plotit=FALSE)$keep
  #  temp<-temp[flag,]
  #}
  x <- temp[,1]
  y <- temp[,2]
  pts <- x
  pts <- as.matrix(pts)
  mat <- matrix(NA,nrow=nboot,ncol=nrow(pts))
  vals <- NA
  for(it in 1:nboot){
    idat <- sample(c(1:length(y)), replace = TRUE)
    xx <- temp[idat,1]
    yy <- temp[idat,2]
    mat[it,] <- runhat(xx,yy,pts=pts,est=est,fr=fr)
  }
  rmd <- apply(mat,2,mean,na.rm=TRUE)
  return(rmd)
}
