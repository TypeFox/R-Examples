dpearsonVII <- function(x,df,location,scale,params,log=FALSE) {
  if (!missing(params)) { df <- params[[1]]; location <- params[[2]]
                          scale <- params[[3]] }
#  stopifnot((df>0)&&(scale!=0))
  if (log) {
    dt((x-location)/scale,df=df,log=TRUE)-log(abs(scale))
  } else {
    1/abs(scale)*dt((x-location)/scale,df=df)
  }  
}

ppearsonVII <- function(q,df,location,scale,params,lower.tail=TRUE,log.p=FALSE) {
  if (!missing(params)) { df <- params[[1]]; location <- params[[2]]
                          scale <- params[[3]] }
#  stopifnot((df>0)&&(scale!=0))
  pt((q-location)/scale,df=df,lower.tail=xor(scale<0,lower.tail),log.p=log.p)
}

qpearsonVII <- function(p,df,location,scale,params,lower.tail=TRUE,log.p=FALSE) {
  if (!missing(params)) { df <- params[[1]]; location <- params[[2]]
                          scale <- params[[3]] }
#  stopifnot((df>0)&&(scale!=0))
  scale*qt(p,df=df,lower.tail=xor(scale<0,lower.tail),log.p=log.p)+location
}

rpearsonVII <- function(n,df,location,scale,params) {
  if (!missing(params)) { df <- params[[1]]; location <- params[[2]]
                          scale <- params[[3]] }
#  stopifnot((df>0)&&(scale!=0))
  scale*rt(n,df=df)+location
}
