pickMeinshausen <- function(p, PM, select = seq_along(p), alpha=0.05, silent=FALSE) {
  
  if (ncol(PM)!=length(p) & nrow(PM)!=length(p)) stop("invalid permutation matrix")
  else if (ncol(PM)==length(p)) PM<-t(PM)
  nperms <- ncol(PM)
  
  sP<-apply(PM,2,sort)
  sPmd=apply(sP,1,sort)
  quantile <- 1
  
  if (alpha==0) too.high <- FALSE else too.high <- TRUE
  while(too.high){
    count=sum(apply(sP,2,function(t) any(t< sPmd[quantile,]))) 
    if (count>(alpha*nperms)) {
      too.high <- FALSE }
    else {
      quantile <- quantile+1
    }
    
  }
  if (quantile==1) crit <- rep(-1,length(p)) 
  else crit<-sPmd[quantile-1,] 
  
  nrej<-length(select)
  sp <- sort(p[select])
  
  lag<-sapply(1:length(sp), function(x) {sum(crit[1:x] > sp[x] )})
  clag <- cummax(lag)
  out<-max(lag)
  
  if (!silent) {
    cat(nrej, " hypotheses selected. At confidence level ", 1-alpha, ":\n", sep="")
    cat("False null-hypotheses >= ", out, "; ", sep="")
    cat("True null-hypotheses <= ", nrej-out, ".\n", sep="")
    invisible(nrej-out)
  } else
    
    out
}


curveMeinshausen <- function(p, PM, select = seq_along(p), order, alpha=0.05, plot = TRUE) {
  
  selected <- !missing(select) || missing(order)
  ordered <- !missing(order)
  
  if (ordered & selected) 
    stop("please provide either select or order, but not both")
  
  if (ncol(PM)!=length(p) & nrow(PM)!=length(p)) stop("invalid permutation matrix")
  else if (ncol(PM)==length(p)) PM<-t(PM)
  nperms <- ncol(PM)
  
  sP<-apply(PM,2,sort)
  sPmd=apply(sP,1,sort)
  quantile <- 1
  
  if (alpha==0) too.high <- FALSE else too.high <- TRUE
  while(too.high){
    count=sum(apply(sP,2,function(t) any(t< sPmd[quantile,]))) 
    if (count>(alpha*nperms)) {
      too.high <- FALSE }
    else {
      quantile <- quantile+1
    }
    
  }
  if (quantile==1) crit <- rep(-1,length(p)) 
  else crit<-sPmd[quantile-1,]   
  
  if (selected) {
    nrej<-length(select)
    sp <- sort(p[select])
    lag<-sapply(1:length(sp), function(x) {sum(crit[1:x] > sp[x] )})
    res <- cummax(lag)
  } else
  {
    nrej<-length(order)
    res <- numeric(length(order))
    for (i in 1:length(order)){
      sp <- sort( p[order[1:i]] )
      lag<-sapply(1:i, function(x) {sum(crit[1:x] > sp[x] )})
      res[i]<-cummax(lag)[i]
    }
  }
  
  if (plot) {
    false <- c(0, res)
    xs <- 1:length(false)-.5
    tots <- 0:length(res)
    plot(xs, tots, type="S", xlab="number of hypotheses", ylab="number of false null-hypotheses", lty=2)
    lines(xs, false, type="S")
    legend("topleft", c(paste("false null-hypotheses (", 100*(1-alpha), "% conf.)", sep=""),"others"), lty=1:2)
    invisible(res)
  } else
    res
}
