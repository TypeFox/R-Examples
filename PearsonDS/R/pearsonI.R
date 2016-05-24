dpearsonI <- function(x,a,b,location,scale,params,log=FALSE) {
  if (!missing(params)) { a <- params[[1]]; b <- params[[2]]
                          location <- params[[3]]; scale <- params[[4]] }
#  stopifnot((a>0)&&(b>0)&&(scale!=0))
  if (log) {
    dbeta((x-location)/scale,a,b,log=TRUE)-log(abs(scale))
  } else {
    1/abs(scale)*dbeta((x-location)/scale,a,b)
  }  
}

ppearsonI <- function(q,a,b,location,scale,params,lower.tail=TRUE,log.p=FALSE) {
  if (!missing(params)) { a <- params[[1]]; b <- params[[2]]
                          location <- params[[3]]; scale <- params[[4]] }
#  stopifnot((a>0)&&(b>0)&&(scale!=0))
  pbeta((q-location)/scale,a,b,lower.tail=xor(scale<0,lower.tail),log.p=log.p)
}

qpearsonI <- function(p,a,b,location,scale,params,lower.tail=TRUE,log.p=FALSE) {
  if (!missing(params)) { a <- params[[1]]; b <- params[[2]]
                          location <- params[[3]]; scale <- params[[4]] }
#  stopifnot((a>0)&&(b>0)&&(scale!=0))
  scale*qbeta(p,a,b,lower.tail=xor(scale<0,lower.tail),log.p=log.p)+location
}

rpearsonI <- function(n,a,b,location,scale,params) {
  if (!missing(params)) { a <- params[[1]]; b <- params[[2]]
                          location <- params[[3]]; scale <- params[[4]] }
#  stopifnot((a>0)&&(b>0)&&(scale!=0))
  scale*rbeta(n,a,b)+location
}
