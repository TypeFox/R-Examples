## Potential improvement: insert special methods for symmetric beta distribution 
## especially: Inverting the symmetrical beta distribution (->Pierre L'Ecuyer)

dpearsonII <- function(x,a,location,scale,params,log=FALSE) {
  if (!missing(params)) { a <- params[[1]]; 
                          location <- params[[2]]; scale <- params[[3]] }
#  stopifnot((a>0)&&(scale!=0))
  if (log) {
    dbeta((x-location)/scale,a,a,log=TRUE)-log(abs(scale))
  } else {
    1/abs(scale)*dbeta((x-location)/scale,a,a)
  }  
}

ppearsonII <- function(q,a,location,scale,params,lower.tail=TRUE,log.p=FALSE) {
  if (!missing(params)) { a <- params[[1]]; 
                          location <- params[[2]]; scale <- params[[3]] }
#  stopifnot((a>0)&&(scale!=0))
  pbeta((q-location)/scale,a,a,lower.tail=xor(scale<0,lower.tail),log.p=log.p)
}

qpearsonII <- function(p,a,location,scale,params,lower.tail=TRUE,log.p=FALSE) {
  if (!missing(params)) { a <- params[[1]]; 
                          location <- params[[2]]; scale <- params[[3]] }
#  stopifnot((a>0)&&(scale!=0))
  scale*qbeta(p,a,a,lower.tail=xor(scale<0,lower.tail),log.p=log.p)+location
}

rpearsonII <- function(n,a,location,scale,params) {
  if (!missing(params)) { a <- params[[1]]; 
                          location <- params[[2]]; scale <- params[[3]] }
#  stopifnot((a>0)&&(scale!=0))
  scale*rbeta(n,a,a)+location
}
