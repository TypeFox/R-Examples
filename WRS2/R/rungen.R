rungen<-function(x, y, fr = 1, est = "mom"){
  #
  # running  interval smoother that can  be used  with any measure
  # of location or scale. By default, an M-estimator is used.
  #
  # LP=TRUE, the plot is further smoothed via lows
  #
  # fr controls amount of smoothing
  eout <- FALSE
  xout <- FALSE
  pyhat <- TRUE
  
  est <- match.arg(est, c("mom", "onestep", "median"), several.ok = FALSE)
  est <- get(est)
  m <- cbind(x,y)
  m <- elimna(m)
  #if(eout && xout)stop("Not allowed to have eout=xout=T")
  #if(eout){
  #  flag<-outfun(m,plotit=FALSE)$keep
  #  m<-m[flag,]
  #}
  #if(xout){
  #  flag<-outfun(m[,1])$keep
  #  m<-m[flag,]
  #}
  x <- m[,1]
  y <- m[,2]
  rmd <- c(1:length(x))
  for(i in 1:length(x)) rmd[i] <- est(y[near(x,x[i],fr)])
  return(rmd)  
}
