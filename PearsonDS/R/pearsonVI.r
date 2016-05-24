dpearsonVI <- function(x,a,b,location,scale,params,log=FALSE) {
  if (!missing(params)) { a <- params[[1]]; b <- params[[2]]
                          location <- params[[3]]; scale <- params[[4]] }
#  stopifnot((a>0)&&(b>0)&&(scale!=0))
  nscale <- scale*a/b
  if (log) {
    df((x-location)/nscale,2*a,2*b,log=TRUE)-log(abs(nscale))
  } else {
    1/abs(nscale)*df((x-location)/nscale,2*a,2*b)
  }  
}

ppearsonVI <- function(q,a,b,location,scale,params,lower.tail=TRUE,log.p=FALSE) {
  if (!missing(params)) { a <- params[[1]]; b <- params[[2]]
                          location <- params[[3]]; scale <- params[[4]] }
#  stopifnot((a>0)&&(b>0)&&(scale!=0))
  nscale <- scale*a/b
  pf((q-location)/nscale,2*a,2*b,lower.tail=xor(scale<0,lower.tail),log.p=log.p)
}

qpearsonVI <- function(p,a,b,location,scale,params,lower.tail=TRUE,log.p=FALSE) {
  if (!missing(params)) { a <- params[[1]]; b <- params[[2]]
                          location <- params[[3]]; scale <- params[[4]] }
#  stopifnot((a>0)&&(b>0)&&(scale!=0))
  nscale <- scale*a/b
  nscale*qf(p,2*a,2*b,lower.tail=xor(scale<0,lower.tail),log.p=log.p)+location
}

rpearsonVI <- function(n,a,b,location,scale,params) {
  if (!missing(params)) { a <- params[[1]]; b <- params[[2]]
                          location <- params[[3]]; scale <- params[[4]] }
#  stopifnot((a>0)&&(b>0)&&(scale!=0))
  nscale <- scale*a/b
  nscale*rf(n,2*a,2*b)+location
}
