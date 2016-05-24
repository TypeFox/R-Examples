dpearson0 <- function(x,mean,sd,params,log=FALSE) {
  if (!missing(params)) { mean <- params[[1]]; sd <- params[[2]] }
#  stopifnot(sd>0)
  dnorm(x,mean=mean,sd=sd,log=log)
}

ppearson0 <- function(q,mean,sd,params,lower.tail=TRUE,log.p=FALSE) {
  if (!missing(params)) { mean <- params[[1]]; sd <- params[[2]] }
#  stopifnot(sd>0)
  pnorm(q,mean=mean,sd=sd,lower.tail=lower.tail,log.p=log.p)
}

qpearson0 <- function(p,mean,sd,params,lower.tail=TRUE,log.p=FALSE) {
  if (!missing(params)) { mean <- params[[1]]; sd <- params[[2]] }
#  stopifnot(sd>0)
  qnorm(p,mean=mean,sd=sd,lower.tail=lower.tail,log.p=log.p)
}

rpearson0 <- function(n,mean,sd,params) {
  if (!missing(params)) { mean <- params[[1]]; sd <- params[[2]] }
#  stopifnot(sd>0)
  rnorm(n,mean=mean,sd=sd)
}
