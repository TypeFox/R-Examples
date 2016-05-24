MVNorm <- function(loc=c(0,0), scale = diag(length(loc))){
   dim0 <- length(loc)
   param <- new("MVNormParameter", loc = loc, scale = as.matrix(scale),
                 name = gettext("parameter of multivariate normal distribution"))
   img0 <- new("EuclideanSpace", dimension = round(dim0,digits=0))

   sigma <- scale %*% t(scale)

   rfun <- function(n){}
   body(rfun) <- substitute({t(rmvnorm(n, mean = loc0, sigma = sigma0))},
                         list(loc0 = loc, sigma0 = sigma))
   dfun <- function(x, log = FALSE){}
   body(dfun) <- substitute({ dmvnorm(x, mean = loc0, sigma = sigma0, log = log)},
                         list(loc0 = loc, sigma0 = sigma))
   pfun <- function(lower=-Inf, upper=Inf){}
   body(pfun) <- substitute( {
                  if(is.matrix(lower)){
                    nr <- nrow(lower)
                    if(nr!=length(loc0))
                       stop("Number of rows must equal dimension of the distribution.")
                    nc <- ncol(lower)
                    if(is.matrix(upper))
                       if( nrow(upper)!=nr || ncol(upper)!=nc)
                           stop("Mismatch of argument dimensions of 'lower' and 'upper'.")
                  }else{
                    nc <- 1
                    nr <- length(loc0)
                    if(is.matrix(upper)){
                       nr <- nrow(upper)
                       if(nr!=length(loc0))
                          stop("Number of rows must equal dimension of the distribution.")
                       nc <- ncol(upper)
                    }
                  }
                    lower <- matrix(rep(lower, length.out = nr*nc),
                                    ncol = nc, nrow = nr)
                    upper <- matrix(rep(upper, length.out = nr*nc),
                                    ncol = nc, nrow = nr)

                    return(sapply(1:nc, function(i)
                            pmvnorm(lower = lower[,i], upper = upper[,i],
                                    mean = loc0, sigma = sigma0)))},
                         list(loc0 = loc, sigma0 = sigma))
                         
   qfun <- function(p, interval = c(-10, 10), tail = c("lower.tail",
                         "upper.tail", "both.tails")){}
   body(qfun) <- substitute( {qmvnorm(p = p, interval = interval, tail = tail,
                              mean = loc0, sigma = sigma0)},
                         list(loc0 = loc, sigma0 = sigma))

   new("MVNormDistribution",
              radDistr =sqrt(Chisq(df=dim0)),
              param = param , p = pfun, q = qfun, r = rfun, d = dfun,
              img = img0, .withSim = FALSE, .withArith = FALSE,
              .logExact = TRUE, .lowerExact = TRUE,
              Symmetry = EllipticalSymmetry(loc))
   }

## MVNormParameter
setMethod("mean", "MVNormParameter",
           function(x) object@loc)
setMethod("sigma", "MVNormParameter",
           function(object) object@scale%*%t(object@scale))


## MVNormDistribution
setMethod("sigma", "MVNormDistribution",
           function(object) object@param@scale%*%t(object@param@scale))

setMethod("mean", "MVNormDistribution",
           function(x) object@param@loc)
