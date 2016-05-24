SNorm <- function(mean=0, sd=1, xi=1.5){
  if(!length(mean)==1L) stop("Parameter 'mean' must be of length 1.")
  if(!length(sd)==1L) stop("Parameter 'sd' must be of length 1.")
  if(!length(xi)==1L) stop("Parameter 'xi' must be of length 1.")
  if(xi<0) stop("Parameter 'xi' must be positive.")
  if(sd<0) stop("Parameter 'sd' must be positive.")
  m0 <- mean; sd0 <- sd; xi0 <- xi
  rfn <- substitute({rsnorm(n,mean=m0s,sd=sd0s,xi=xi0s)},
                     list(m0s=m0,sd0s=sd0,xi0s=xi0))
  dfn <- substitute({d0 <- dsnorm(x,mean=m0s,sd=sd0s,xi=xi0s)
                     return(if(log) log(d0) else d0)},
                     list(m0s=m0,sd0s=sd0,xi0s=xi0))
  pfn <- substitute({ p00 <- psnorm(q,mean=m0s,sd=sd0s,xi=xi0s)
                      p0  <- if(lower.tail) p00 else 1-p00
                      return(if(log.p) log(p0) else p0)},
                     list(m0s=m0,sd0s=sd0,xi0s=xi0))
  qfn <- substitute({p00 <- if(log.p) exp(p) else p
                          p0 <- if(lower.tail) p00 else 1-p00
                          return(qsnorm(p0,mean=m0s,sd=sd0s,xi=xi0s))},
                     list(m0s=m0,sd0s=sd0,xi0s=xi0))
  rF <- function(n){}; body(rF) <- rfn
  dF <- function(x, log=FALSE){}; body(dF) <- dfn
  pF <- function(q, lower.tail=TRUE, log.p=FALSE){}; body(pF) <- pfn
  qF <- function(p, lower.tail=TRUE, log.p=FALSE){}; body(qF) <- qfn
  new("SNorm",r = rF, d = dF, p = pF, q = qF,
       param = new("SNormParameter", mean = m0,sd = sd0,xi = xi0))
}
STd <- function(mean=0, sd=1, nu=5){
  if(!length(mean)==1L) stop("Parameter 'mean' must be of length 1.")
  if(!length(sd)==1L) stop("Parameter 'sd' must be of length 1.")
  if(!length(nu)==1L) stop("Parameter 'nu' must be of length 1.")
  if(sd<0) stop("Parameter 'sd' must be positive.")
  if(nu<=2L) stop("Parameter 'nu' must be larger than 2.")
  sqrt((nu-2)/nu)*sd*Td(ncp=0,df=nu)+mean
}
SSTd <- function(mean=0, sd=1, nu=5, xi=1.5){
  if(!length(mean)==1L) stop("Parameter 'mean' must be of length 1.")
  if(!length(sd)==1L) stop("Parameter 'sd' must be of length 1.")
  if(!length(nu)==1L) stop("Parameter 'nu' must be of length 1.")
  if(!length(xi)==1L) stop("Parameter 'xi' must be of length 1.")
  if(xi<0) stop("Parameter 'xi' must be positive.")
  if(sd<0) stop("Parameter 'sd' must be positive.")
  if(nu<=2L) stop("Parameter 'nu' must be larger than 2.")
  m0 <- mean; sd0 <- sd; nu0 <- nu; xi0 <- xi
  rfn <- substitute({rsstd(n,mean=m0s,sd=sd0s,nu=nu0s, xi=xi0s)},
                     list(m0s=m0,sd0s=sd0,xi0s=xi0,nu0s=nu0))
  dfn <- substitute({d0 <- dsstd(x,mean=m0s,sd=sd0s,xi=xi0s, nu=nu0s)
                     return(if(log) log(d0) else d0)},
                     list(m0s=m0,sd0s=sd0,xi0s=xi0,nu0s=nu0))
  pfn <- substitute({ p00 <- psstd(q,mean=m0s,sd=sd0s,xi=xi0s, nu=nu0s)
                      p0  <- if(lower.tail) p00 else 1-p00
                      return(if(log.p) log(p0) else p0)},
                     list(m0s=m0,sd0s=sd0,xi0s=xi0,nu0s=nu0))
  qfn <- substitute({p00 <- if(log.p) exp(p) else p
                          p0 <- if(lower.tail) p00 else 1-p00
                          return(qsstd(p0,mean=m0s,sd=sd0s,xi=xi0s, nu=nu0s))},
                     list(m0s=m0,sd0s=sd0,xi0s=xi0,nu0s=nu0))
  rF <- function(n){}; body(rF) <- rfn
  dF <- function(x, log=FALSE){}; body(dF) <- dfn
  pF <- function(q, lower.tail=TRUE, log.p=FALSE){}; body(pF) <- pfn
  qF <- function(p, lower.tail=TRUE, log.p=FALSE){}; body(qF) <- qfn
  new("SSTd", r = rF, d = dF, p = pF, q = qF,
       param = new("SSTdParameter", mean=m0,sd=sd0,nu=nu0,xi=xi0)
       )
}

## Access methods
setMethod("xi", signature(object = "SNormParameter"), function(object) object@xi)
setMethod("xi", signature(object = "SSTdParameter"), function(object) object@xi)
setMethod("nu", signature(object = "SSTdParameter"), function(object) object@nu)
setMethod("mean", signature(x = "SSTdParameter"), function(x, ...) x@mean)
setMethod("sd", signature(x = "SSTdParameter"), function(x, ...) x@sd)
## wrapped access methods
setMethod("mean", signature(x = "SNorm"), function(x, ...) mean(param(x)))
setMethod("mean", signature(x = "SSTd"), function(x, ...) mean(param(x)))
setMethod("sd", signature(x = "SNorm"), function(x) sd(param(x)))
setMethod("sd", signature(x = "SSTd"), function(x) sd(param(x)))
setMethod("xi", signature(object = "SNorm"), function(object) object@param@xi)
setMethod("xi", signature(object = "SSTd"), function(object) object@param@xi)
setMethod("nu", signature(object = "SSTd"), function(object) object@param@nu)
## wrapped replace methods
setMethod("mean<-", "SNorm",
           function(object, value) SNorm(mean=value, sd=sd(object), xi=xi(object)))
setMethod("mean<-", "SSTd",
           function(object, value) SSTd(mean=value, sd=sd(object), xi=xi(object)))
setMethod("sd<-", "SNorm",
           function(object, value) SNorm(mean=mean(object), sd=value, xi=xi(object)))
setMethod("sd<-", "SSTd",
           function(object, value) SSTd(mean=mean(object), sd=value, xi=xi(object), nu =nu(object)))
setMethod("nu<-", "SSTd",
           function(object, value) SSTd(mean=mean(object), sd=sd(object), nu=value, xi=xi(object)))
setMethod("xi<-", "SNorm",
           function(object, value) SNorm(mean=mean(object), sd=sd(object), xi=value))
setMethod("xi<-", "SSTd",
           function(object, value) SSTd(mean=mean(object), sd=sd(object), xi=value, nu =nu(object)))
