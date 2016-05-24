# Function from package 'CompQuadForm' v.1.4.1 (c) 2013 P. Lafaye de Micheaux

davies <- function(q,lambda,h = rep(1,length(lambda)),delta = rep(0,length(lambda)),sigma=0,lim=10000,acc=0.0001) {

  r <- length(lambda)
  if (length(h) != r) stop("lambda and h should have the same length!")
  if (length(delta) != r) stop("lambda and delta should have the same length!")
  
  out <- .C("qfc",lambdas=as.double(lambda),noncentral=as.double(delta),df=as.integer(h),r=as.integer(r),sigma=as.double(sigma),q=as.double(q),lim=as.integer(lim),acc=as.double(acc),trace=as.double(rep(0,7)),ifault=as.integer(0),res=as.double(0),PACKAGE="FREGAT")

  out$res <- 1 - out$res
  
  return(list(trace=out$trace,ifault=out$ifault,Qq=out$res))
  
}

