#### as to whether to use Generating functions or to use initialize methods:
#### http://tolstoy.newcastle.edu.au/R/e2/devel/07/01/1976.html
                     
################################################################################
## SPACES
################################################################################

setMethod("initialize", "Reals",
          function(.Object) {
            .Object@dimension <-  1
            .Object@name <- gettext("Real Space")
            .Object
          })


setMethod("initialize", "Naturals",
          function(.Object) {
            .Object@dimension <-  1
            .Object@name <- gettext("Grid of Naturals")
            .Object
          })


################################################################################
## PARAMETERS
################################################################################

setMethod("initialize", "GeomParameter",
          function(.Object, prob = .5) {
            .Deprecated(new = "new(\"NbinomParameter\"(size = 1, prob, name)",
                        package = "distr", 
                        msg = gettext(
"Class 'GeomParameter' is no longer needed and will be replaced by \nclass 'NbinomParameter' soon."                        
                        ))
            .Object@prob <- prob
            .Object@name <- gettext("Parameter of a Geometric distribution")
            .Object
          })
################################################################################
## DISTRIBUTIONS
################################################################################

## Class: UnivariateDistribution
###produces difficulties in coercing...:
#
#setMethod("initialize", "UnivariateDistribution",
#          function(.Object, r = NULL, d = NULL, p = NULL, q = NULL, 
#                    param = NULL, img = new("Reals"),
#                    .withSim = FALSE, .withArith = FALSE) {
#            if(is.null(r)) {
#              stop("You have at least to give the slot r.")
#              return(invisible())}
#            ### Attention: no checking!!!
#            .Object@img <- img
#            .Object@param <- param
#            .Object@d <- d
#            .Object@p <- p
#            .Object@q <- q
#           .Object@r <- r
#            .Object@.withSim <- .withSim
#            .Object@.withArith <- .withArith
#            .Object })

## class AbscontDistribution
setMethod("initialize", "AbscontDistribution",
          function(.Object, r = NULL, d = NULL, p = NULL, q = NULL, 
                   gaps = NULL, param = NULL, img = new("Reals"),
                   .withSim = FALSE, .withArith = FALSE,
                   .lowerExact = FALSE, .logExact = FALSE,
                   low1 = NULL, up1 = NULL, low = -Inf, up =Inf,
                   Symmetry = NoSymmetry()
                   ) {
            ## don't use this if the call is new("AbscontDistribution")
            LL <- length(sys.calls())
            if(sys.calls()[[LL-3]] == "new(toDef)")
               {return(.Object)}
            if(sys.calls()[[LL-3]] == "new(\"AbscontDistribution\")")
               {return(.Object)}
            
            if(is.null(r))
               warning("you have to specify slot r at least")
                          
            ## TOBEDONE Errorkanal
            
            dpq.approx <- 0
            
            dfun <- d
            pfun <- p
            qfun <- q
            
            if(is.null(d)) {
              .withSim <- TRUE
              dpq <- RtoDPQ(r)
              dpq.approx <- 1
              dfun <- dpq$dfun}
            
            if(is.null(p)) {
              .withSim <- TRUE
              if(dpq.approx == 0) {dpq <- RtoDPQ(r)}
              dpq.approx <- 1
              pfun <- dpq$pfun}
            
            if(is.null(q)) {
               ## quantile function
               rN <- NULL
               if(is.null(up1)) up1 <- max(rN <- r(10^getdistrOption("RtoDPQ.e")))
               if(is.null(low1)) {
                 low1 <- if(is.null(rN)) min(r(10^getdistrOption("RtoDPQ.e")))
                         else min(rN)}
                         
               h <- (up1-low1)/getdistrOption("DefaultNrFFTGridPointsExponent")
               x <-   seq(from = low1, to = up1, by = h)

               px.l <- pfun(x + 0.5*h)
               px.u <- pfun(x + 0.5*h, lower.tail = FALSE)
            
               qfun <- .makeQNew(x + 0.5*h, px.l, px.u, FALSE, low, up)
             }
            
            .Object@img <- img
            .Object@param <- param
            .Object@d <- dfun
            .Object@p <- pfun
            .Object@q <- qfun
            .Object@r <- r
            .Object@gaps <- gaps
            .Object@.withSim <- .withSim
            .Object@.withArith <- .withArith
            .Object@.logExact <- .logExact
            .Object@.lowerExact <- .lowerExact
            .Object@Symmetry <- Symmetry
            .Object })

## class AffLinAbscontDistribution
setMethod("initialize", "AffLinAbscontDistribution",
          function(.Object, r = NULL, d = NULL, p = NULL, q = NULL, gaps = NULL,
                   a = 1, b = 0, X0 = Norm(), param = NULL, img = new("Reals"),
                   .withSim = FALSE, .withArith = FALSE,
                   .lowerExact = FALSE, .logExact = FALSE,
                   Symmetry = NoSymmetry()) {
  X <- new("AbscontDistribution", r = r, d = d, p = p, q = q, gaps = gaps, 
           param = param, img = img, .withSim = .withSim, 
           .withArith = .withArith)
  .Object@gaps  <- X@gaps 
  .Object@img <- X@img
  .Object@param <- X@param
  .Object@a <- a
  .Object@b <- b
  .Object@X0 <- X0
  .Object@d <- X@d
  .Object@p <- X@p
  .Object@q <- X@q
  .Object@r <- X@r
  .Object@.withSim <- .withSim
  .Object@.withArith <- .withArith
  .Object@Symmetry <- Symmetry
  .Object
})

## Class: DiscreteDistribution
setMethod("initialize", "DiscreteDistribution",
          function(.Object, r = NULL, d = NULL, p = NULL, q = NULL, 
                    support = NULL, param = NULL, img = new("Reals"), 
                    .withSim = FALSE, .withArith = FALSE,
                   .lowerExact = FALSE, .logExact = FALSE,
                   Symmetry = NoSymmetry()) {

            ## don't use this if the call is new("DiscreteDistribution")
            LL <- length(sys.calls())
            if(sys.calls()[[LL-3]] == "new(toDef)")
               {return(.Object)}
            if(sys.calls()[[LL-3]] == "new(\"DiscreteDistribution\")")
               {return(.Object)}
            
            if(is.null(r))
               warning("you have to specify slot r at least")
              
            if(is.null(support)) 
               .Object@support <- as.numeric(names(table(r(10^6))))
            
            else .Object@support <- support

#           len = length(support)
#
#            if(len > 1){
#              if(min(diff(support)) < getdistrOption("DistrResolution"))
#                 stop("grid too narrow --> change DistrResolution")
#            }

            dpq.approx <- 0

            dfun <- d
            pfun <- p
            qfun <- q

            if(is.null(d)) {
              .withSim <- TRUE
              dpq <- RtoDPQ.d(r)
              dpq.approx <- 1
              dfun <- dpq$dfun
            }

            if(is.null(p)) {
              .withSim <- TRUE
              if(dpq.approx==0) dpq <- RtoDPQ.d(r)
              dpq.approx <- 1
              pfun <- dpq$pfun
            }

            if(is.null(q)) {
              .withSim <- TRUE
              if(dpq.approx==0) dpq <- RtoDPQ.d(r)
              qfun <- dpq$qfun
            }

            .Object@img <- img
            .Object@param <- param
            .Object@d <- dfun
            .Object@p <- pfun
            .Object@q <- qfun
            .Object@r <- r
            .Object@.withSim <- .withSim
            .Object@.withArith <- .withArith
            .Object@.lowerExact <- .lowerExact
            .Object@.logExact <- .logExact
            .Object@Symmetry <- Symmetry
            .Object
          })

## Class: AffLinDiscreteDistribution
setMethod("initialize", "AffLinDiscreteDistribution",
          function(.Object, r = NULL, d = NULL, p = NULL, q = NULL, 
                   support = NULL, a = 1, b = 0, X0 = Binom(), param = NULL, 
                   img = new("Reals"), .withSim = FALSE, .withArith = FALSE,
                   .lowerExact = FALSE, .logExact = FALSE,
                   Symmetry = NoSymmetry()) {
   ## don't use this if the call is new("DiscreteDistribution")
   LL <- length(sys.calls())
   if(sys.calls()[[LL-3]] == "new(\"AffLinDiscreteDistribution\")")
        X <- new("DiscreteDistribution")
   else X <- new("DiscreteDistribution", r = r, d = d, p = p, q = q, support = support, 
             param = param, img = img, .withSim = .withSim, 
            .withArith = .withArith)
  .Object@support  <- X@support 
  .Object@img <- X@img
  .Object@param <- X@param
  .Object@a <- a
  .Object@b <- b
  .Object@X0 <- X0
  .Object@d <- X@d
  .Object@p <- X@p
  .Object@q <- X@q
  .Object@r <- X@r
  .Object@.withSim <- .withSim
  .Object@.withArith <- .withArith
  .Object@.lowerExact <- .lowerExact
  .Object@.logExact <- .logExact
  .Object@Symmetry <- Symmetry
  .Object
})

## Class: LatticeDistribution
setMethod("initialize", "LatticeDistribution",
          function(.Object, r = NULL, d = NULL, p = NULL, q = NULL, 
                    support = NULL, lattice = NULL, param = NULL, 
                    img = new("Reals"), .withSim = FALSE, .withArith = FALSE,
                   .lowerExact = FALSE, .logExact = FALSE,
                   Symmetry = NoSymmetry()) {


             LL <- length(sys.calls())
             if(sys.calls()[[LL-3]] == "new(\"LatticeDistribution\")")
             D <- new("DiscreteDistribution")
             else
             D <- new("DiscreteDistribution", r = r, d = d, p = p, 
                       q = q, support = support, param = param, img = img, 
                     .withSim = .withSim, .withArith = .withArith)

            
             OS  <- D@support 

             #if(is.null(lattice))  
             #  {  if(! .is.vector.lattice(OS))
             #         stop("Support as given/generated is not a lattice.")
             #     .Object@lattice <- .make.lattice.es.vector(OS)
             #}else{
                  .Object@lattice <- if(is.null(lattice )) 
                          new("Lattice") else lattice
             #}


            .Object@support <- OS
            .Object@img <- D@img
            .Object@param <- D@param
            .Object@d <- D@d
            .Object@p <- D@p
            .Object@q <- D@q
            .Object@r <- D@r
            .Object@.withSim <- .withSim
            .Object@.withArith <- .withArith
            .Object@.lowerExact <- .lowerExact
            .Object@.logExact <- .logExact
            .Object@Symmetry <- Symmetry
            .Object
          })

## Class: AffLinLatticeDistribution
setMethod("initialize", "AffLinLatticeDistribution",
          function(.Object, r = NULL, d = NULL, p = NULL, q = NULL, 
                   support = NULL, lattice = NULL, a = 1, b = 0, X0 = Binom(), 
                   param = NULL, img = new("Reals"), .withSim = FALSE, 
                   .withArith = FALSE, .lowerExact = FALSE, .logExact = FALSE,
                   Symmetry = NoSymmetry()) {

   LL <- length(sys.calls())
   if(sys.calls()[[LL-3]] == "new(\"AffLinLatticeDistribution\")")
        X <- new("LatticeDistribution")
   else X <- new("LatticeDistribution", r = r, d = d, p = p, q = q, 
                  support = support, lattice = lattice, param = param, 
                  img = img, .withSim = .withSim, 
                 .withArith = .withArith)

  .Object@support  <- X@support 
  .Object@lattice <-  X@lattice 
  .Object@img <- X@img
  .Object@param <- X@param
  .Object@a <- a
  .Object@b <- b
  .Object@X0 <- X0
  .Object@d <- X@d
  .Object@p <- X@p
  .Object@q <- X@q
  .Object@r <- X@r
  .Object@.withSim <- .withSim
  .Object@.withArith <- .withArith
  .Object@.lowerExact <- .lowerExact
  .Object@.logExact <- .logExact
  .Object@Symmetry <- Symmetry
  .Object
})

######### particular discrete distributions

### Class: Dirac distribution
setMethod("initialize", "Dirac",
          function(.Object, location = 0, .withArith = FALSE) {
            .Object@img <- new("Reals")
            .Object@param <- new("DiracParameter", location = location)
            .Object@r <- function(n){}
            .Object@d <- function(x, log = FALSE){} 
            .Object@p <- function(q, lower.tail = TRUE, log.p = FALSE){} 
            .Object@q <- function(p, lower.tail = TRUE, log.p = FALSE){}
            body(.Object@r) <- substitute({ rep(locationSub, n)},
                                            list(locationSub = location)
                                          )
            body(.Object@d) <- substitute(
                           { y <- rep(locationSub, length(x))
                             d0 <- mapply(function(x,y) 
                                          as.numeric(isTRUE(all.equal(x,y))),
                                          x = x, y = y)
                             if (log) d0 <- log(d0)
                             return(d0)
                           }, list(locationSub = location)
                                          )
            body(.Object@p) <- substitute(
                           {p0 <-as.numeric(q + .Machine$double.eps^.5 >= 
                                       locationSub)
                            if (!lower.tail) p0 <- 1-p0
                            if (log.p) p0 <- log(p0)
                            return(p0)
                            },
                            list(locationSub = location)
                                          )
            body(.Object@q) <- substitute( 
                { if (log.p) p <- exp(p)
                  if(any((p < 0)|(p > 1))) 
                     warning("q Method of class Dirac produced NaNs.")
                  q0 <- ifelse((p < 0)|(p > 1), NaN, locationSub) 
                  return(q0)
                },
                           list(locationSub = location)
                                          )
            .Object@support <- location
            .Object@lattice <- new("Lattice", pivot = location, width = 1, 
                                    Length = 1)
            .Object@.withArith <- .withArith
            .Object
          })

## Class: binomial distribution
setMethod("initialize", "Binom",
          function(.Object, size = 1, prob = 0.5, .withArith = FALSE) {
            .Object@img <- new("Naturals")
            .Object@param <- new("BinomParameter", size = size, prob = prob)
            .Object@support <- 0:size
            .Object@r <- function(n){}
            .Object@d <- function(x, log = FALSE){}
            .Object@p <- function(q, lower.tail = TRUE, log.p = FALSE){} 
            .Object@q <- function(p, lower.tail = TRUE, log.p = FALSE){} 
            body(.Object@r) <- substitute(
                           { rbinom(n, size = sizeSub, prob = probSub) },
                             list(sizeSub = size, probSub = prob)
                                         )
            body(.Object@d) <- substitute(
                           { dbinom(x, size = sizeSub, prob = probSub, 
                                    log = log) },
                             list(sizeSub = size, probSub = prob)
                                         )
            body(.Object@p) <- substitute(
                           { pbinom(q, size = sizeSub, prob = probSub, 
                                    lower.tail = lower.tail, log.p = log.p) },
                             list(sizeSub = size, probSub = prob)
                                         )
            body(.Object@q) <- substitute(
                           { qbinom(p, size = sizeSub, prob = probSub, 
                                    lower.tail = lower.tail, log.p = log.p) },
                             list(sizeSub = size, probSub = prob)
                                         )
            .Object@support = 0:size
            .Object@lattice = new("Lattice", pivot = 0, width = 1, 
                                   Length = size+1)
            .Object@.withArith <- .withArith
            .Object
          })

## Class: hypergeometric distribution
setMethod("initialize", "Hyper",
          function(.Object, m = 1, n = 1, k = 1, .withArith = FALSE) {
            .Object@img <- new("Naturals")
            .Object@param <- new("HyperParameter", m = m, n = n, k = k)
            .Object@support <- 0:k
            .Object@r <- function(nn){}
            .Object@d <- function(x, log = FALSE){}
            .Object@p <- function(q, lower.tail = TRUE, log.p = FALSE){} 
            .Object@q <- function(p, lower.tail = TRUE, log.p = FALSE){} 
            body(.Object@r) <- substitute(
                               { rhyper(nn, m = mSub, n = nSub, k = kSub) },
                                 list(mSub = m, nSub = n, kSub = k)
                                         )
            body(.Object@d) <- substitute(
                               { dhyper(x, m = mSub, n = nSub, k = kSub, 
                                        log = log) },
                                 list(mSub = m, nSub = n, kSub = k)
                                          )
            body(.Object@p) <- substitute(
                               { phyper(q, m = mSub, n = nSub, k = kSub, 
                                        lower.tail = lower.tail, log.p = log.p) 
                                        },
                                 list(mSub = m, nSub = n, kSub = k)
                                          )
            body(.Object@q) <- substitute(
                               { qhyper(p, m = mSub, n = nSub, k = kSub, 
                                        lower.tail = lower.tail, log.p = log.p) 
                                        },
                                 list(mSub = m, nSub = n, kSub = k)
                                          )
            .Object@support <-  seq(from = 0, to = min(k,m), by = 1)
            .Object@lattice <-  new("Lattice", pivot = 0, width = 1,
                                     Length = min(k,m)+1 )
            .Object@.withArith <- .withArith
            .Object
          })

## Class: Poisson distribution 
setMethod("initialize", "Pois",
          function(.Object, lambda = 1, .withArith=FALSE) {
            .Object@img <- new("Naturals")
            .Object@param <- new("PoisParameter", lambda = lambda)
            .Object@r <- function(n){}
            .Object@d <- function(x, log = FALSE){}
            .Object@p <- function(q, lower.tail = TRUE, log.p = FALSE){} 
            .Object@q <- function(p, lower.tail = TRUE, log.p = FALSE){} 
            body(.Object@r) <- substitute({ rpois(n, lambda = lambdaSub) },
                                            list(lambdaSub = lambda)
                                          )
            body(.Object@d) <- substitute({ dpois(x, lambda = lambdaSub, 
                                                  log = log) },
                                            list(lambdaSub = lambda)
                                          )
            body(.Object@p) <- substitute({ ppois(q, lambda = lambdaSub, 
                                                  lower.tail = lower.tail, 
                                                  log.p = log.p) },
                                            list(lambdaSub = lambda)
                                          )
            body(.Object@q) <- substitute({ qpois(p, lambda = lambdaSub, 
                                                  lower.tail = lower.tail, 
                                                  log.p = log.p) },
                                            list(lambdaSub = lambda)
                                          )
            .Object@support <- seq(from = 0, by = 1, to = 
                                   qpois(getdistrOption("TruncQuantile"),
                                         lambda = lambda, lower.tail = FALSE) 
                                         + 2
                                   )
            .Object@lattice <- new("Lattice", pivot = 0, width = 1, 
                                    Length = Inf)
            .Object@.withArith <- .withArith
            .Object
          })

## Class: negative binomial distribution
setMethod("initialize", "Nbinom",
          function(.Object, size = 1, prob = 0.5, .withArith = FALSE) {
            .Object@img <- new("Naturals")
            .Object@param <- new("NbinomParameter", size = size, prob = prob)
            .Object@r <- function(n){}
            .Object@d <- function(x, log = FALSE){}
            .Object@p <- function(q, lower.tail = TRUE, log.p = FALSE){} 
            .Object@q <- function(p, lower.tail = TRUE, log.p = FALSE){} 
            body(.Object@r) <- substitute(
                           { rnbinom(n, size = sizeSub, prob = probSub) },
                             list(sizeSub = size, probSub = prob)
                                         )
            body(.Object@d) <- substitute(
                           { dnbinom(x, size = sizeSub, prob = probSub, 
                                     log = log) },
                             list(sizeSub = size, probSub = prob)
                                         )
            body(.Object@p) <- substitute(
                           { pnbinom(q, size = sizeSub, prob = probSub, 
                                     lower.tail = lower.tail, log.p = log.p) },
                             list(sizeSub = size, probSub = prob)
                                         )
            body(.Object@q) <- substitute(
                           { qnbinom(p, size = sizeSub, prob = probSub, 
                                     lower.tail = lower.tail, log.p = log.p) },
                             list(sizeSub = size, probSub = prob)
                                         )
            .Object@.withArith <- .withArith
            .Object@support <-  seq(from = 0, by = 1, 
                                    to = qnbinom( size = size, prob = prob,
                                         getdistrOption("TruncQuantile"),
                                         lower.tail = FALSE)
                                    )
            .Object@lattice <-  new("Lattice", pivot = 0, width = 1, 
                                     Length = Inf)
            .Object
          })

## Class: geometric distribution
setMethod("initialize", "Geom",
          function(.Object, prob = 0.5, .withArith = FALSE) {
            .Object@img <- new("Naturals")
            .Object@param <- new("NbinomParameter", name = 
                             gettext("Parameter of a Geometric distribution"),
                             prob = prob)
            .Object@support <- 0:qgeom(getdistrOption("TruncQuantile"), 
                                       prob = prob, lower.tail = FALSE)
            .Object@r <- function(n){}
            .Object@d <- function(x, log = FALSE){}
            .Object@p <- function(q, lower.tail = TRUE, log.p = FALSE){} 
            .Object@q <- function(p, lower.tail = TRUE, log.p = FALSE){} 
            body(.Object@r) <- substitute({ rgeom(n, prob = probSub) },
                                          list(probSub = prob))
            body(.Object@d) <- substitute({ dgeom(x, prob = probSub, 
                                                  log = log) },
                                          list(probSub = prob))
            body(.Object@p) <- substitute({ pgeom(q, prob = probSub, 
                                                  lower.tail = lower.tail, 
                                                  log.p = log.p) },
                                          list(probSub = prob))
            body(.Object@q) <- substitute({ qgeom(p, prob = probSub, 
                                                  lower.tail = lower.tail, 
                                                  log.p = log.p) },
                                          list(probSub = prob))
            .Object@.withArith <- .withArith
            .Object
          })


## --- particular absolutely continuous distributions


## Class: uniform distribution
setMethod("initialize", "Unif",
          function(.Object, Min = 0, Max = 1, .withArith = FALSE) {
            .Object@img <- new("Reals")
            .Object@param <- new("UnifParameter", Min = Min, Max = Max)
            .Object@r <- function(n){}
            .Object@d <- function(x, log = FALSE){}
            .Object@p <- function(q, lower.tail = TRUE, log.p = FALSE){} 
            .Object@q <- function(p, lower.tail = TRUE, log.p = FALSE){} 
            body(.Object@r) <- substitute(
                           { runif(n, min = MinSub, max = MaxSub) },
                             list(MinSub = Min, MaxSub = Max)
                                          )
            body(.Object@d) <- substitute(
                           { dunif(x, min = MinSub, max = MaxSub, log = log) },
                             list(MinSub = Min, MaxSub = Max)
                                          )            
            body(.Object@p) <- substitute(
                           { punif(q, min = MinSub, max = MaxSub, 
                                   lower.tail = lower.tail, log.p = log.p) },
                             list(MinSub = Min, MaxSub = Max)
                                          )        
            body(.Object@q) <- substitute(
                           { qunif(p, min = MinSub, max = MaxSub, 
                                   lower.tail = lower.tail, log.p = log.p) },
                             list(MinSub = Min, MaxSub = Max)
                                          )                    
            .Object@.withArith <- .withArith
            .Object@Symmetry <- SphericalSymmetry(Min+Max/2)
            .Object
          })


## Class: normal distribution
setMethod("initialize", "Norm",
          function(.Object, mean = 0, sd = 1, .withArith = FALSE) {
            .Object@img <- new("Reals")
            .Object@param <- new("UniNormParameter", mean = mean, sd = sd)
            .Object@r <- function(n){}
            .Object@d <- function(x, log = FALSE){}
            .Object@p <- function(q, lower.tail = TRUE, log.p = FALSE){} 
            .Object@q <- function(p, lower.tail = TRUE, log.p = FALSE){} 
            body(.Object@r) <- substitute(
                           { rnorm(n, mean = meanSub, sd = sdSub) },
                             list(meanSub = mean, sdSub = sd)
                                          )
            body(.Object@d) <- substitute(
                           { dnorm(x, mean = meanSub, sd = sdSub, log = log) },
                             list(meanSub = mean, sdSub = sd)
                                          )
            body(.Object@p) <- substitute(
                           { pnorm(q, mean = meanSub, sd = sdSub, 
                                   lower.tail = lower.tail, log.p = log.p) },
                             list(meanSub = mean, sdSub = sd)
                                          )
            body(.Object@q) <- substitute(
                           { qnorm(p, mean = meanSub, sd = sdSub, 
                                   lower.tail = lower.tail, log.p = log.p) },
                             list(meanSub = mean, sdSub = sd)
                                          )
            .Object@.withArith <- .withArith
            .Object@Symmetry <- SphericalSymmetry(mean)
            .Object
          })

## Class: lognormal distribution
setMethod("initialize", "Lnorm",
          function(.Object, meanlog = 0, sdlog = 1, .withArith = FALSE) {
            .Object@img <- new("Reals")
            .Object@param <- new("LnormParameter", meanlog = meanlog, 
                                  sdlog = sdlog)
            .Object@r <- function(n){}
            .Object@d <- function(x, log = FALSE){}
            .Object@p <- function(q, lower.tail = TRUE, log.p = FALSE){} 
            .Object@q <- function(p, lower.tail = TRUE, log.p = FALSE){} 
            body(.Object@r) <- substitute(
                           { rlnorm(n, meanlog = meanlogSub, 
                                    sdlog = sdlogSub) },
                             list(meanlogSub = meanlog, sdlogSub = sdlog)
                                          )
            body(.Object@d) <- substitute(
                           { dlnorm(x, meanlog = meanlogSub, 
                                    sdlog = sdlogSub, log = log) },
                             list(meanlogSub = meanlog, sdlogSub = sdlog)
                                          )
            body(.Object@p) <- substitute(
                           { plnorm(q, meanlog = meanlogSub, sdlog = sdlogSub, 
                                    lower.tail = lower.tail, log.p = log.p) },
                             list(meanlogSub = meanlog, sdlogSub = sdlog)
                                          )
            body(.Object@q) <- substitute(
                           { qlnorm(p, meanlog = meanlogSub, sdlog = sdlogSub, 
                                    lower.tail = lower.tail, log.p = log.p) },
                             list(meanlogSub = meanlog, sdlogSub = sdlog)
                                          )
            .Object@.withArith <- .withArith
            .Object
          })

## Class: CauchyDistribution
setMethod("initialize", "Cauchy",
          function(.Object, location = 0, scale = 1) {
            .Object@img <- new("Reals")
            .Object@param <- new("CauchyParameter", location = location, 
                                  scale = scale)
            .Object@r <- function(n){}
            .Object@d <- function(x, log = FALSE){}
            .Object@p <- function(q, lower.tail = TRUE, log.p = FALSE){} 
            .Object@q <- function(p, lower.tail = TRUE, log.p = FALSE){} 
            body(.Object@r) <- substitute(
                           { rcauchy(n, location = locationSub, 
                                     scale = scaleSub) },
                             list(locationSub = location, 
                                               scaleSub = scale)
                                          )
            body(.Object@d) <- substitute(
                           { dcauchy(x, location = locationSub, 
                                     scale = scaleSub, log = log) }, 
                             list(locationSub = location, scaleSub = scale)
                                          )
            body(.Object@p) <- substitute(
                           { pcauchy(q, location = locationSub, 
                                     scale = scaleSub, lower.tail = lower.tail, 
                                     log.p = log.p) },
                             list(locationSub = location, scaleSub = scale)
                                          )
            body(.Object@q) <- substitute(
                           { qcauchy(p, location = locationSub, 
                                     scale = scaleSub, lower.tail = lower.tail, 
                                     log.p = log.p) },
                             list(locationSub = location, scaleSub = scale)
                                          )
            .Object@.withSim   <- FALSE
            .Object@.withArith <- FALSE
            .Object@Symmetry <- SphericalSymmetry(location)
            .Object
          })

## Class: F distribution
setMethod("initialize", "Fd",
          function(.Object, df1 = 1, df2 = 1, ncp = 0) {
            .Object@img <- new("Reals")
            .Object@param <- new("FParameter", df1 = df1, df2 = df2, ncp = ncp)
            .Object@r <- function(n){}            
            .Object@d <- function(x, log = FALSE){}
            .Object@p <- function(q, lower.tail = TRUE, log.p = FALSE){} 
            .Object@q <- function(p, lower.tail = TRUE, log.p = FALSE){} 
            
            #### will probably change.... (when df for ncp!=0 available...)
            if((isTRUE(all.equal(ncp,0)))||getRversion()>'2.4.0')
               {df.0<- function(x, df1 = df1, df2 = df2, ncp = ncp, log = FALSE)
                       {stats::df(x = x, df1 = df1 , df2 = df2, ncp = ncp, 
                                  log = log)} 
               }
            else  
               {## for R < 2.4.0  df with ncp != 0:
                ### later perhaps with sfsmisc:
                  TQ <- getdistrOption("TruncQuantile")/2
                  xz <- qf(TQ, df1 = df1, df2 = df2, ncp = ncp, 
                           lower.tail = FALSE)
                  pfun <- function(x){pf(x, df1 = df1, df2 = df2, ncp = ncp)}
                  dfun <- .P2D(p=pfun, ql = 0, qu = xz)
                # by means of simulations
                # rfun <- function(n){rf(n, df1=df1, df2=df2, ncp=ncp)}
                # dfun <-R2D(rfun, nsim = 10^getdistrOption("RtoDPQ.e"), 
                #           n = getdistrOption("DefaultNrGridPoints"))
                df.0 <- function(x, df1 = df1, df2 = df2, ncp = ncp, 
                                 log = FALSE) {dfun(x)}
               }                  
            
            body(.Object@r) <- substitute(                    
                           { rf(n, df1 = df1Sub, df2 = df2Sub, ncp = ncpSub) },
                             list(df1Sub = df1, df2Sub = df2, ncpSub = ncp)                                         
                                          )
            body(.Object@d) <- substitute(
                           { df.0(x, df1 = df1Sub, df2 = df2Sub, ncp = ncpSub, 
                                  log = log)},
                             list(df1Sub = df1, df2Sub = df2, ncpSub = ncp)
                                           )
            body(.Object@p) <- substitute(
                           { pf(q, df1 = df1Sub, df2 = df2Sub, ncp = ncpSub, 
                                lower.tail = lower.tail, log.p = log.p) },
                             list(df1Sub = df1, df2Sub = df2, ncpSub = ncp)
                                          )
            body(.Object@q) <- substitute(
                           { qf(p, df1 = df1Sub, df2 = df2Sub, ncp = ncpSub, 
                                lower.tail = lower.tail, log.p = log.p) },
                             list(df1Sub = df1, df2Sub = df2, ncpSub = ncp)
                                          )
            .Object@.withArith <- FALSE
            .Object
          })

## Class: Student distribution
setMethod("initialize", "Td",
          function(.Object, df = 1, ncp = 0) {
            .Object@img <- new("Reals")
            .Object@param <- new("TParameter", df = df, ncp = ncp) 
            .Object@r <- function(n){}
            .Object@d <- function(x, log = FALSE){}
            .Object@p <- function(q, lower.tail = TRUE, log.p = FALSE){} 
            .Object@q <- function(p, lower.tail = TRUE, log.p = FALSE){} 
            body(.Object@r) <- substitute({ rt(n, df = dfSub, ncp = ncpSub) }, 
                                          list(dfSub = df, ncpSub = ncp)
                                          )
            body(.Object@d) <- substitute(
                                       { dt(x, df = dfSub, ncp = ncpSub, 
                                            log = log) },
                                         list(dfSub = df, ncpSub = ncp)
                                          )
            body(.Object@p) <- substitute(
                                       { pt(q, df = dfSub, ncp = ncpSub, 
                                            lower.tail = lower.tail, 
                                            log.p = log.p) },
                                         list(dfSub = df, ncpSub = ncp)
                                          )
            body(.Object@q) <- substitute(
                                       { qt(p, df = dfSub, ncp = ncpSub, 
                                            lower.tail = lower.tail, 
                                            log.p = log.p) },
                                         list(dfSub = df, ncpSub = ncp)
                                          )
            .Object@.withArith <- FALSE
            .Object@Symmetry <- SphericalSymmetry(0)
            .Object
          })

## Class: Chi squared distribution
setMethod("initialize", "Chisq",
          function(.Object, df = 1, ncp = 0, .withArith = FALSE) {
            .Object@img <- new("Reals")
            .Object@param <- new("ChisqParameter", df = df, ncp = ncp)
            .Object@r <- function(n){}
            .Object@d <- function(x, log = FALSE){}
            .Object@p <- function(q, lower.tail = TRUE, log.p = FALSE){} 
            .Object@q <- function(p, lower.tail = TRUE, log.p = FALSE){} 
            body(.Object@r) <- substitute(
                           { rchisq(n, df = dfSub, ncp = ncpSub) },
                             list(dfSub = df, ncpSub = ncp)
                                          )
            body(.Object@d) <- substitute(
                           { dchisq(x, df = dfSub, ncp = ncpSub, log = log) },
                             list(dfSub = df, ncpSub = ncp)
                                          )
            body(.Object@p) <- substitute(
                           { pchisq(q, df = dfSub, ncp = ncpSub, 
                                    lower.tail = lower.tail, log.p = log.p) },
                             list(dfSub = df, ncpSub = ncp)
                                          )
            body(.Object@q) <- substitute(
                           { qchisq(p, df = dfSub, ncp = ncpSub, 
                                    lower.tail = lower.tail, log.p = log.p) },
                             list(dfSub = df, ncpSub = ncp)
                                          )
            .Object@.withSim   <- FALSE
            .Object@.withArith <- .withArith
            .Object
          })

## Class: exponential distribution
setMethod("initialize", "Exp",
          function(.Object, rate = 1, .withArith = FALSE) {
            .Object@img <- new("Reals")
            .Object@param <- new("ExpParameter", rate = rate)
            .Object@r <- function(n){}
            .Object@d <- function(x, log = FALSE){}
            .Object@p <- function(q, lower.tail = TRUE, log.p = FALSE){} 
            .Object@q <- function(p, lower.tail = TRUE, log.p = FALSE){} 
            body(.Object@r) <- substitute(
                           { rexp(n, rate = rateSub) },
                             list(rateSub = rate)
                                          )
            body(.Object@d) <- substitute(
                           { dexp(x, rate = rateSub, log = log) },
                             list(rateSub = rate)
                                          )
            body(.Object@p) <- substitute(
                           { pexp(q, rate = rateSub, 
                                  lower.tail = lower.tail, log.p = log.p) },
                             list(rateSub = rate)
                                          )
            body(.Object@q) <- substitute(
                           { qexp(p, rate = rateSub, 
                                  lower.tail = lower.tail, log.p = log.p) },
                             list(rateSub = rate)
                                          )
            .Object@.withSim   <- FALSE
            .Object@.withArith <- .withArith
            .Object
          })

## Class: Laplace or Double Exponential distribution
setMethod("initialize", "DExp",
          function(.Object, rate = 1, .withArith = FALSE) {
            .Object@img <- new("Reals")
            .Object@param <- new("ExpParameter", rate = rate)
            .Object@r <- function(n){}
            .Object@d <- function(x, log = FALSE){}
            .Object@p <- function(q, lower.tail = TRUE, log.p = FALSE){} 
            .Object@q <- function(p, lower.tail = TRUE, log.p = FALSE){} 
            body(.Object@r) <- substitute(
                           { (2 * rbinom(n, size = 1, prob = 0.5) -1 ) * 
                              rexp(n, rate = rateSub) 
                           }, list(rateSub = rate)
                                          )
            body(.Object@d) <- substitute( 
                            { d0 <-  dexp(abs(x), rate = rateSub, log = log) 
                              d0 <- if (log) d0-log(2) else d0 <- d0 / 2
                              return(d0) },
                              list(rateSub = rate)
                                          )
            body(.Object@p) <- substitute(
                           { if (!lower.tail) q <- -q
                             p0 <- ifelse( q <= 0, 
                                           0.5 * pexp(-q, rate = rateSub,
                                                    lower.tail = FALSE),
                                           0.5 + 0.5*pexp( q, rate = rateSub)
                                           )
                             if (log.p)  p0 <- log(p0)       
                             return(p0)
                           }, list(rateSub = rate)
                                         )
            body(.Object@q) <- substitute(
                           {  if (log.p) p <- exp(p)
                              if (!lower.tail) p <- 1-p
                              ifelse( p <= 0.25,          
                                  -qexp(2*p, rate = rateSub, lower.tail =FALSE),
                                  ifelse( p <= 0.5,
                                      -qexp(1-2*p, rate = rateSub),
                                      ifelse( p <= 0.75   ,
                                          qexp(2*p - 1, rate = rateSub),
                                          qexp(2*(1-p), rate = rateSub, 
                                               lower.tail = FALSE) 
                                            ) 
                                         ) 
                                     )
                           }, list(rateSub = rate)
                                          )
            .Object@.withSim   <- FALSE
            .Object@.withArith <- .withArith
            .Object@Symmetry <- SphericalSymmetry(0)
            .Object
          })

## Class: gamma distribution
setMethod("initialize", "Gammad",
          function(.Object, shape = 1, scale = 1, .withArith = FALSE) {
            .Object@img <- new("Reals")
            .Object@param <- new("GammaParameter", shape = shape, scale = scale)
            .Object@r <- function(n){}
            .Object@d <- function(x, log = FALSE){}
            .Object@p <- function(q, lower.tail = TRUE, log.p = FALSE){} 
            .Object@q <- function(p, lower.tail = TRUE, log.p = FALSE){} 
            body(.Object@r) <- substitute(
                           { rgamma(n, shape = shapeSub, scale = scaleSub) },
                             list(shapeSub = shape, scaleSub = scale)
                                          )
            body(.Object@d) <- substitute(
                           { dgamma(x, shape = shapeSub, scale = scaleSub, 
                                    log = log) },
                             list(shapeSub = shape, scaleSub = scale)
                                          )
            body(.Object@p) <- substitute(
                           { pgamma(q, shape = shapeSub, scale = scaleSub, 
                                    lower.tail = lower.tail, log.p = log.p) },
                             list(shapeSub = shape, scaleSub = scale)
                                          )
            body(.Object@q) <- substitute(
                           { qgamma(p, shape = shapeSub, scale = scaleSub, 
                                    lower.tail = lower.tail, log.p = log.p) },
                             list(shapeSub = shape, scaleSub = scale)
                                          )
            .Object@.withSim   <- FALSE
            .Object@.withArith <- .withArith
            .Object
          })


## Class: BetaDistribution
setMethod("initialize", "Beta",
          function(.Object, shape1 = 1, shape2 = 1, ncp = 0) {
            .Object@img <- new("Reals")
            .Object@param <- new("BetaParameter", shape1 = shape1, 
                                  shape2 = shape2, ncp = ncp)
            .Object@r <- function(n){}
            .Object@d <- function(x, log = FALSE){}
            .Object@p <- function(q, lower.tail = TRUE, log.p = FALSE){} 
            .Object@q <- function(p, lower.tail = TRUE, log.p = FALSE){} 
            body(.Object@r) <- substitute(
                           { rbeta(n, shape1 = shape1Sub, shape2 = shape2Sub, 
                                   ncp = ncpSub) },
                             list(shape1Sub = shape1, shape2Sub = shape2, 
                                  ncpSub = ncp)
                                          )
            body(.Object@d) <- substitute(
                           { dbeta(x, shape1 = shape1Sub, shape2 = shape2Sub, 
                                   ncp = ncpSub, log = log) },
                             list(shape1Sub = shape1, shape2Sub = shape2, 
                                  ncpSub = ncp)
                                          )
            body(.Object@p) <- substitute(
                           { pbeta(q, shape1 = shape1Sub, shape2 = shape2Sub, 
                                   ncp = ncpSub, lower.tail = lower.tail, 
                                   log.p = log.p) },
                             list(shape1Sub = shape1, shape2Sub = shape2, 
                                  ncpSub = ncp)
                                          )
            body(.Object@q) <- substitute(
                           { qbeta(p, shape1 = shape1Sub, shape2 = shape2Sub, 
                                   ncp = ncpSub, lower.tail = lower.tail, 
                                   log.p = log.p) },
                             list(shape1Sub = shape1, shape2Sub = shape2, 
                                  ncpSub = ncp)
                                          )
            .Object@.withSim   <- FALSE
            .Object@.withArith <- FALSE
            .Object
          })

## Class: logistic distribution
setMethod("initialize", "Logis",
          function(.Object, location = 0, scale = 1) {
            .Object@img <- new("Reals")
            .Object@param <- new("LogisParameter", location = location, 
                                  scale = scale)
            .Object@r <- function(n){}
            .Object@d <- function(x, log = FALSE){}
            .Object@p <- function(q, lower.tail = TRUE, log.p = FALSE){} 
            .Object@q <- function(p, lower.tail = TRUE, log.p = FALSE){} 
            body(.Object@r) <- substitute(
                           { rlogis(n, location = locationSub, 
                                    scale = scaleSub) },
                             list(locationSub = location, scaleSub = scale)
                                         )
            body(.Object@d) <- substitute(
                           { dlogis(x, location = locationSub, scale = scaleSub, 
                                    log = log) },
                             list(locationSub = location, scaleSub = scale)
                                         )
            body(.Object@p) <- substitute(
                           { plogis(q, location = locationSub, scale = scaleSub, 
                                    lower.tail = lower.tail, log.p = log.p) },
                             list(locationSub = location, scaleSub = scale)
                                         )
            body(.Object@q) <- substitute(
                           { qlogis(p, location = locationSub, scale = scaleSub, 
                                    lower.tail = lower.tail, log.p = log.p) },
                             list(locationSub = location, scaleSub = scale)
                                         )
            .Object@.withArith <- FALSE
            .Object
          })

## Class: Weibull distribution
setMethod("initialize", "Weibull",
          function(.Object, shape = 1, scale = 1) {
            .Object@img <- new("Reals")
            .Object@param <- new("WeibullParameter", 
                                  shape = shape, scale = scale
                                  )
            .Object@r <- function(n){}
            .Object@d <- function(x, log = FALSE){}
            .Object@p <- function(q, lower.tail = TRUE, log.p = FALSE){} 
            .Object@q <- function(p, lower.tail = TRUE, log.p = FALSE){} 
            body(.Object@r) <- substitute(
                           { rweibull(n, shape = shapeSub, scale = scaleSub) },
                             list(shapeSub = shape, scaleSub = scale)
                                         )
            body(.Object@d) <- substitute(
                           { dweibull(x, shape = shapeSub, scale = scaleSub, 
                                      log = log) },
                             list(shapeSub = shape, scaleSub = scale)
                                         )
            body(.Object@p) <- substitute(
                           { pweibull(q, shape = shapeSub, scale = scaleSub, 
                                      lower.tail = lower.tail, log.p = log.p) },
                             list(shapeSub = shape, scaleSub = scale)
                                         )
            body(.Object@q) <- substitute(
                           { qweibull(p, shape = shapeSub, scale = scaleSub, 
                                      lower.tail = lower.tail, log.p = log.p) },
                             list(shapeSub = shape, scaleSub = scale)
                                         )
            .Object@.withArith <- FALSE
            .Object
          })

## Class: Arcsine distribution
setMethod("initialize", "Arcsine",
          function(.Object, .withArith = FALSE) {
            .Object@img <- new("Reals")
            .Object@r <- function(n){sin((runif(n)-.5)*pi)}
            .Object@d <- function(x, log = FALSE){ 
                              x0 <- (abs(x)<1-.Machine$double.eps)
                              x1 <- x^2*x0
                              d <-  x0/sqrt(1-x1)/pi
                              d[.isEqual(abs(x),1)] <- Inf
                              if(log) d<- log(d)
                              return(d)}
            .Object@p <- function(q, lower.tail = TRUE, log.p = FALSE){ 
                              if(!lower.tail) q<- -q
                              q <- pmin(pmax(q,-1),1)
                              p <- asin(q)/pi+1/2
                              if(log.p) p <- log(p)
                              return(p)} 
            .Object@q <- function(p, lower.tail = TRUE, log.p = FALSE){ 
                              if(log.p) p <- exp(p)
                              p1 <- p
                              p1[p<0|p>1] <- 0.5
                              if(!lower.tail) p1 <- 1-p1
                              q <- sin( (p1-1/2)*pi)
                              q[p<0|p>1] <- NA
                              q[.isEqual(p,0)] <- -1
                              q[.isEqual(p,1)] <-  1
                              return(q)}                      
            .Object@.withSim   <- FALSE
            .Object@.withArith <- .withArith
            .Object@Symmetry <- SphericalSymmetry(0)
            .Object
          })
