MVt <- function(loc = c(0,0), scale = diag(length(loc)), df = 1, ncp = 0){

   dim0 <- length(loc)
   img0 <- new("EuclideanSpace", dimension = round(dim0,digits=0))
   param <- new("MVtParameter", loc=loc, scale = as.matrix(scale),
                 df = df, ncp = ncp)
   sigma <- scale %*% t(scale)

   rfun <- function(n){}
   body(rfun) <- substitute(t(rmvnorm(n, sigma = sigma0)/
                                (rchisq(n, df = df0, ncp = ncp0)/df0)^.5
                              )+loc0,
                         list(loc0 = loc, sigma0 = sigma, df0=df, ncp0=ncp))
   dfun <- function(x, log = FALSE){}
   body(dfun) <- substitute( dmvt(t(x - loc0), sigma = sigma0, df = df0,
                                  delta = ncp0, log = log),
                         list(loc0 = loc, sigma0 = sigma, df0=df, ncp0=ncp))
   pfun <- function(lower=-Inf, upper=Inf){}
   body(pfun) <- substitute( pmvt(lower=t(lower-loc0), upper=t(upper-loc0),
                                  sigma = sigma0, df = df0, delta = ncp0),
                         list(loc0 = loc, sigma0 = sigma, df0=df, ncp0=ncp))
   qfun <- function(p, interval = c(-10, 10), tail = c("lower.tail",
                         "upper.tail", "both.tails")){}
   body(qfun) <- substitute({q0 <-  qmvt(p = p, interval = interval, tail = tail,
                              sigma = sigma0, df = df0, delta = ncp0)
                              q0$quantile <- q0$quantile+loc0
                              q0},
                         list(loc0 = loc, sigma0 = sigma, df0=df, ncp0=ncp))

   dradfun <- function(x, log = FALSE){
                x0 <- x
                x0[x<0] <- 1
                lg <- lgamma((dim0+df)/2)-lgamma(df/2)+(dim0-1)*log(x) -
                      lgamma(dim0/2) +log(2) - dim0/2*log(df)-
                      (dim0+df)/2*log(1+x^2/df)
                lg[x<0] <- -Inf
                return(if(log) lg else exp(lg))}

   radDistr <- AbscontDistribution(d = dradfun)

   new("MVtDistribution",
        radDistr = radDistr,
        param = param ,p =pfun, q=qfun, r=rfun, d=dfun,
        img = img0, .withSim = FALSE, .withArith = FALSE,
        .logExact = TRUE, .lowerExact = TRUE,
        Symmetry = EllipticalSymmetry(loc))
   }

## MVtParameter
setMethod("sigma", "MVtParameter",
           function(object) object@scale%*%t(object@scale))
setMethod("df", "MVtParameter", function(x) x@df)
setMethod("ncp", "MVtParameter", function(object) object@ncp)

## MVtDistribution
setMethod("sigma", "MVtDistribution",
           function(object) object@param@scale%*%t(object@param@scale))

setMethod("df", "MVtDistribution", function(x) x@param@df)

setMethod("ncp", "MVtDistribution", function(object) object@param@ncp)

