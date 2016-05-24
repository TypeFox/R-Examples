dpearsonIII <- function(x,shape,location,scale,params,log=FALSE) {
  if (!missing(params)) { shape <- params[[1]]; 
                          location <- params[[2]]; scale <- params[[3]] }
#  stopifnot((shape>0)&&(scale!=0))
  gscale <- abs(scale); ssgn <- sign(scale)
  dgamma(ssgn*(x-location),shape=shape,scale=gscale,log=log)
}

ppearsonIII <- function(q,shape,location,scale,params,lower.tail=TRUE,
                        log.p=FALSE) {
  if (!missing(params)) { shape <- params[[1]]; 
                          location <- params[[2]]; scale <- params[[3]] }
#  stopifnot((shape>0)&&(scale!=0))
  gscale <- abs(scale); ssgn <- sign(scale)
  pgamma(ssgn*(q-location),shape=shape,scale=gscale,
         lower.tail=xor(scale<0,lower.tail),log.p=log.p)
}

qpearsonIII <- function(p,shape,location,scale,params,lower.tail=TRUE,
                        log.p=FALSE) {
  if (!missing(params)) { shape <- params[[1]]; 
                          location <- params[[2]]; scale <- params[[3]] }
#  stopifnot((shape>0)&&(scale!=0))
  gscale <- abs(scale); ssgn <- sign(scale)
  ssgn*qgamma(p,shape=shape,scale=gscale,
              lower.tail=xor(scale<0,lower.tail),log.p=log.p)+location
}

rpearsonIII <- function(n,shape,location,scale,params) {
  if (!missing(params)) { shape <- params[[1]]; 
                          location <- params[[2]]; scale <- params[[3]] }
#  stopifnot((shape>0)&&(scale!=0))
  gscale <- abs(scale); ssgn <- sign(scale)
  ssgn*rgamma(n,shape=shape,scale=gscale)+location
}
