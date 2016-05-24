dpearsonV <- function(x,shape,location,scale,params,log=FALSE) {
  if (!missing(params)) { shape <- params[[1]]; 
                          location <- params[[2]]; scale <- params[[3]] }
#  stopifnot((shape>0)&&(scale!=0))
  gscale <- abs(scale); ssgn <- sign(scale)
  if (log) {
    tmp <- dgamma(1/pmax(ssgn*(x-location),0),shape=shape,scale=gscale,log=TRUE)-
             2*log(pmax(ssgn*(x-location),0))
    tmp[is.na(tmp)] <- -Inf
    tmp         
  } else {
    ifelse(x==location,0,dgamma(1/pmax(ssgn*(x-location),0),shape=shape,scale=gscale)/((x-location)^2))
  }  
}

ppearsonV <- function(q,shape,location,scale,params,lower.tail=TRUE,log.p=FALSE) {
  if (!missing(params)) { shape <- params[[1]]; 
                          location <- params[[2]]; scale <- params[[3]] }
#  stopifnot((shape>0)&&(scale!=0))
  gscale <- abs(scale); ssgn <- sign(scale)
  ifelse(q==location,as.numeric(!(scale>0)),
         pgamma(1/pmax(ssgn*(q-location),0),shape=shape,scale=gscale,
           lower.tail=xor(scale>0,lower.tail),log.p=log.p))
}

qpearsonV <- function(p,shape,location,scale,params,lower.tail=TRUE,log.p=FALSE) {
  if (!missing(params)) { shape <- params[[1]]; 
                          location <- params[[2]]; scale <- params[[3]] }
#  stopifnot((shape>0)&&(scale!=0))
  gscale <- abs(scale); ssgn <- sign(scale)
  ssgn/qgamma(p,shape=shape,scale=gscale,lower.tail=xor(scale>0,lower.tail),
              log.p=log.p)+location
}

rpearsonV <- function(n,shape,location,scale,params) {
  if (!missing(params)) { shape <- params[[1]]; 
                          location <- params[[2]]; scale <- params[[3]] }
#  stopifnot((shape>0)&&(scale!=0))
  gscale <- abs(scale); ssgn <- sign(scale)
  ssgn/rgamma(n,shape=shape,scale=gscale)+location
}
