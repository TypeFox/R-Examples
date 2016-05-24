farebrother <- function(q,lambda,h = rep(1,length(lambda)),delta = rep(0,length(lambda)),maxit=100000,eps=10^(-10),mode=1) {


  r <- length(lambda)
  if (length(h) != r) stop("lambda and h should have the same length!")
  if (length(delta) != r) stop("lambda and delta should have the same length!")

#  mode <- 1 # ??

  dnsty <- 0.0
  ifault <- 0
  res <- 0 
  

  out <- .C("ruben",lambda=as.double(lambda),h=as.integer(h),delta=as.double(delta),r=as.integer(r),q=as.double(q),mode=as.double(mode),maxit=as.integer(maxit),eps=as.double(eps),dnsty=as.double(dnsty),ifault=as.integer(ifault),res=as.double(res),PACKAGE="CompQuadForm")

  out$res <- 1 - out$res

  return(out)

}

