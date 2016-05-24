############# preparations ################
# (.onload, .onattach ...)
### left out P.R. 01-09-11
##.onLoad <- function(lib, pkg) { # extended 03-28-06: P.R. 
##    require("methods", character = TRUE, quietly = TRUE)
##}



.onAttach <- function(library, pkg)
{
  unlockBinding(".distroptions", asNamespace("distr"))
  unlockBinding(".distrExInstalled", asNamespace("distr"))

## global variable needed for flat.R
##  unlockBinding(".OkTyp", asNamespace("distr"))
    msga <- gettext(
    "Attention: Arithmetics on distribution objects are understood as "
                   )
    msgb <- gettext(
    "operations on corresponding random variables (r.v.s); see distrARITH().\n"
                   )
    msgc <- gettext(
    "Some functions from package 'stats' are intentionally masked ---see distrMASK().\n"
                   )
    msgd <- gettext(
    "Note that global options are controlled by distroptions() ---c.f. ?\"distroptions\"."
                   )
buildStartupMessage(pkg = "distr", msga, msgb, msgc, msgd, library = library, 
                    packageHelp = TRUE, 
# MANUAL = "http://www.uni-bayreuth.de/departments/math/org/mathe7/DISTR/distr.pdf",
                    VIGNETTE = gettext(
"Package \"distrDoc\" provides a vignette to this package as well as to several extension packages; try vignette(\"distr\")."
                                      )
                   )
  invisible()
} 

################################
##
## Optional..-classes
##
################################
setClassUnion("OptionalMatrix", 
               c("matrix","NULL")
               )
### from Matthias' thesis / ROptEst
## optional numeric
setClassUnion("OptionalNumeric", c("numeric", "NULL"))

################################
##
## utility classes 
##
################################

setClass("Integer", contains ="numeric",
          validity = function(object) all(.isInteger(object)))

################################
##
## space classes 
##
################################

## virtal Class: rSpace
setClass("rSpace", 
          representation = representation(name = "character"), 
          prototype = prototype(name = gettext("a space")), 
          contains = "VIRTUAL"
          )

## Class: EuclideanSpace
setClass("EuclideanSpace", 
          representation = representation(dimension = "numeric"), 
          contains = "rSpace",
          prototype = prototype(dimension = 1, 
                                name = gettext("Euclidean Space")
                                )
         )

## Class: Reals
setClass("Reals",  
          contains = "EuclideanSpace"
          )


## Class: Lattice
setClass("Lattice", 
          representation = representation(pivot = "numeric", width = "numeric", 
                                          Length = "numeric"
### masking not possible here -> Length instead of length
                                          ),
          prototype = prototype(pivot = 0, width = 1, Length = 2, 
                                name = gettext("a default lattice")
                                ),
          contains = "rSpace"
         )

## Class: Naturals
setClass("Naturals", 
          contains = "Reals"
          )

################################
##
## parameter classes
##
################################

setClass("Parameter", 
          representation = representation(name = "character"), 
          prototype = prototype(name = gettext("a parameter"))
          )

setClassUnion("OptionalParameter", 
               c("Parameter","NULL")
               )


## Class: ChisqParameter
setClass("ChisqParameter", 
          representation = representation(df = "numeric", ncp = "numeric"), 
          prototype = prototype(df = 1, ncp = 0, 
                      name = gettext("Parameter of a Chisq distribution")
                      ), 
          contains = "Parameter"
          )

### Class: DiracParameter
setClass("DiracParameter", 
          representation = representation(location = "numeric"), 
          prototype = prototype(location = 0, 
                      name = gettext("Parameter of a Dirac distribution")
                      ), 
          contains = "Parameter"
          )

## Class: ExpParameter
setClass("ExpParameter", 
          representation = representation(rate = "numeric"), 
          prototype = prototype(rate = 1, name = 
                      gettext("Parameter of an Exponential distribution")
                      ),  
          contains = "Parameter"
          )

## Class: GammaParameter
setClass("GammaParameter", 
          representation = representation(shape = "numeric", scale = "numeric"), 
          prototype = prototype(shape = 1, scale = 1, 
                      name = gettext("Parameter of a Gamma distribution")
                      ), 
          contains = "Parameter"
          )

## Class: PoisParameter
setClass("PoisParameter", 
          representation = representation(lambda = "numeric"), 
          prototype = prototype(lambda = 1, 
                      name = gettext("Parameter of a Poisson distribution")
                      ), 
          contains = "Parameter"
          )

## Class: NbinomParameter
setClass("NbinomParameter", 
          representation = representation(size = "numeric", prob = "numeric"), 
          prototype = prototype(size = 1, prob = 0.5, name = 
                      gettext("Parameter of a Negative Binomial distribution")
                      ), 
          contains = "Parameter"
          )

## Class: HyperParameter
setClass("HyperParameter", 
          representation = representation(m = "numeric", n = "numeric", 
                                          k = "numeric"
                                          ), 
          prototype = prototype(m = 1, n = 1, k = 1, name = 
                      gettext("Parameter of a Hypergeometric distribution")
                      ), 
          contains = "Parameter"
          )

## Class: BinomParameter
setClass("BinomParameter", 
          representation = representation(size = "numeric", prob = "numeric"), 
          prototype = prototype(size = 1, prob = 0.5, name = 
                      gettext("Parameter of a Binomial distribution")
                      ), 
          contains = "Parameter"
          )

#-
## no longer needed: this is a negBinom with size 1 no longer 
#-
### !!! deprecated as of version 1.9 !!!
##
## Class: GeomParameter   
setClass("GeomParameter", 
          representation = representation(prob = "numeric"), 
          prototype = prototype(prob = 0.5, name = 
                      gettext("Parameter of a Geometric distribution")
                      ), 
          contains = "Parameter"
          )
### !!! end of deprecated !!! 

## Class: CauchyParameter
setClass("CauchyParameter", 
          representation = representation(location = "numeric", 
                                          scale = "numeric"
                                          ), 
          prototype = prototype(location = 0, scale = 1, name = 
                      gettext("Parameter of a Cauchy distribution")
                      ), 
          contains = "Parameter"
          )

## Class: NormParameter
setClass("NormParameter", 
          representation = representation(mean = "numeric", sd = "vector"), 
          prototype = prototype(mean = 0, sd = 1, name = 
                      gettext("Parameter of a Normal distribution")
                      ), 
          contains = "Parameter"
          )

## Class: UniNormParameter
setClass("UniNormParameter", 
          prototype = prototype(name = 
                      gettext("Parameter of a univariate Normal distribution")
                      ), 
          contains = "NormParameter"
          )

## Class: UnifParameter
setClass("UnifParameter", 
          representation = representation(Min = "numeric", Max = "numeric"), 
          prototype = prototype(Min = 0, Max = 1, name =  
                      gettext("Parameter of a Uniform distribution")
                      ), 
          contains = "Parameter"
          )

## Class: FParameter
setClass("FParameter", 
          representation = representation(df1 = "numeric", df2 = "numeric", 
                                          ncp = "numeric"
                                          ), 
          prototype = prototype(df1 = 1, df2 = 1, ncp = 0, name = 
                      gettext("Parameter of a Fisher-Snedecor/F distribution")
                      ), 
          contains = "Parameter"
          )

## Class: TParameter
setClass("TParameter", 
          representation = representation(df = "numeric", ncp = "numeric"), 
          prototype = prototype(df = 1, ncp = 0, name = 
                      gettext("Parameter of a Student/T distribution")
                      ), 
          contains = "Parameter"
          )

## Class: LNormParameter
setClass("LnormParameter", 
          representation = representation(meanlog = "numeric",
                                          sdlog = "numeric"
                                          ), 
          prototype = prototype(meanlog = 0, meansd = 1, name =  
                      gettext("Parameter of a Log-Normal distribution")
                      ), 
          contains = "Parameter"
          )

## Class: BetaParameter
setClass("BetaParameter", 
          representation = representation(shape1 = "numeric", 
                                          shape2 = "numeric", ncp = "numeric"
                                          ), 
          prototype = prototype(shape1 = 1, shape2 = 1, ncp = 0, name = 
                      gettext("Parameter of a Beta distribution")
                      ), 
          contains = "Parameter"
          )

## Class: LogisParameter
setClass("LogisParameter", 
          representation = representation(location = "numeric", 
                                          scale = "numeric"
                                          ), 
          prototype = prototype(location = 0, scale = 1, name = 
                      gettext("Parameter of a Logistic distribution")
                      ), 
          contains = "Parameter"
          )

## Class: WeibullParameter
setClass("WeibullParameter", 
          representation = representation(shape = "numeric", 
                                          scale = "numeric"
                                          ), 
          prototype = prototype(shape = 1, scale = 1, name = 
                      gettext("Parameter of a Weibull distribution")
                      ), 
          contains = "Parameter"
          )

################################
##
## matrix classes
##
################################
## positive definite, symmetric matrices with finite entries
setClass("PosSemDefSymmMatrix", contains = "matrix",
            prototype = prototype(matrix(1)),
            validity = function(object){
                if(nrow(object) != ncol(object))
                    stop("no square matrix")
                if(any(!is.finite(object)))
                    stop("inifinite or missing values in matrix")
                if(!isTRUE(all.equal(object, t(object), .Machine$double.eps^0.5,
                                     check.attributes = FALSE)))
                    stop("matrix is not symmetric")
                if(!all(eigen(object)$values > -100*.Machine$double.eps))
                   stop("matrix is (numerically) not positive semi - definite")
               return(TRUE)
            })

## positive definite, symmetric matrices with finite entries
setClass("PosDefSymmMatrix", contains = "PosSemDefSymmMatrix",
            validity = function(object){
               if(!all(eigen(object)$values > 100*.Machine$double.eps))
                   stop("matrix is (numerically) not positive definite")
               valid <- getValidity(getClass("PosSemDefSymmMatrix"))
               valid(as(object, "PosSemDefSymmMatrix"))
               return(TRUE)
            })


################################
##
## symmetry classes
##
################################

### from Matthias' thesis / ROptEst / moved from distrMod

## class of symmetries
setClass("Symmetry", representation(type = "character",
                                    SymmCenter = "ANY"),
                     contains = "VIRTUAL")

## symmetry of distributions
setClass("DistributionSymmetry", contains = c("Symmetry", "VIRTUAL"))

## no symmetry
setClass("NoSymmetry", contains = "DistributionSymmetry",
            prototype = prototype(type = "non-symmetric distribution",
                                  SymmCenter = NULL))

## elliptical symmetry
setClass("EllipticalSymmetry", contains = "DistributionSymmetry",
            prototype = prototype(type = "elliptically symmetric distribution",
                                  SymmCenter = numeric(0)))

## spherical symmetry
setClass("SphericalSymmetry", contains = "EllipticalSymmetry",
            prototype = prototype(type = "spherically symmetric distribution",
                                  SymmCenter = numeric(0)))

## list of symmetry types
setClass(Class = "DistrSymmList",
            prototype = prototype(list(new("NoSymmetry"))),
            contains = "list",
            validity = function(object){
                nrvalues <- length(object)
                for(i in 1:nrvalues)
                    if(!is(object[[i]], "DistributionSymmetry"))
                        stop("element ", i, " is no 'DistributionSymmetry'")
                return(TRUE)
            })


################################
##
## distribution classes
##
################################

setClass("Distribution",
          representation = representation(
                      img = "rSpace",
                      param = "OptionalParameter",
                      r = "function",
                      d = "OptionalFunction",
                      p = "OptionalFunction",
                      q = "OptionalFunction", # extended by P.R. 28-03-06
                      .withSim = "logical",   ## 'internal' slots => no
                      .withArith = "logical",  ## accessor/replacement functions
                      .logExact = "logical",
                      .lowerExact = "logical",
                      Symmetry = "DistributionSymmetry"
                      ),
         prototype = prototype(
                     r = function(n){ rnorm(n, mean = 0, sd = 1) },
                     d = function(x, log = FALSE)
                            { dnorm(x, mean = 0, sd = 1, log = log) },
                     p = function(q, lower.tail = TRUE, log.p = FALSE )
                             { pnorm(q, mean = 0, sd = 1,
                                     lower.tail = lower.tail, log.p = log.p) },
                     q = function(p, lower.tail = TRUE, log.p = FALSE )
                             { qnorm(p, mean = 0, sd = 1,
                                     lower.tail = lower.tail, log.p = log.p) },
                     img = new("Reals"),
                     param = NULL,
                     .withArith = FALSE,
                     .withSim = FALSE,
                     .logExact = FALSE,
                     .lowerExact = FALSE,
                     Symmetry = new("NoSymmetry")
                     )
         )

## Class: UnivariateDistribution
setClass("UnivariateDistribution",  
          contains = "Distribution"
          )

### ---- absolutely continuous distributions ---- ###

## AbscontDistribution
setClass("AbscontDistribution", 
          representation = representation(gaps = "OptionalMatrix"),
          prototype = prototype(gaps = NULL),
          contains = "UnivariateDistribution"
          )



## inbetween-Class: ExpOrGammaOrChisq

#not quite virtual ...
setClass("ExpOrGammaOrChisq", 
          contains = c("AbscontDistribution", "VIRTUAL")
          )


## Class: exponential distribution
setClass("Exp",
          prototype = prototype(
                      r = function(n){ rexp(n, rate = 1) },
                      d = function(x, log = FALSE)
                                  { dexp(x, rate = 1, log = log) },
                      p = function(q, lower.tail = TRUE, log.p = FALSE )
                                  { pexp(q, rate = 1, lower.tail = lower.tail,
                                         log.p = log.p) },
                      q = function(q, lower.tail = TRUE, log.p = FALSE )
                                  { qexp(p, rate = 1, lower.tail = lower.tail,
                                         log.p = log.p) },
                      param = new("ExpParameter"),
                     .logExact = TRUE,
                     .lowerExact = TRUE
                      ),
          contains = "ExpOrGammaOrChisq"
          )

## Class: gamma distribution
setClass("Gammad",
          prototype = prototype(
                      r = function(n){ rgamma(n, shape = 1, scale = 1) },
                      d = function(x, log = FALSE){
                              dgamma(x, shape = 1, scale = 1, log = log)
                                          },
                      p = function(q, lower.tail = TRUE, log.p = FALSE ){
                              pgamma(q, shape = 1, scale = 1,
                                     lower.tail = lower.tail, log.p = log.p)
                                          },
                      q = function(p, lower.tail = TRUE, log.p = FALSE ){
                              qgamma(p, shape = 1, scale = 1,
                                     lower.tail = lower.tail, log.p = log.p)
                                          },
                      param = new("GammaParameter"),
                     .logExact = TRUE,
                     .lowerExact = TRUE
                      ),
          contains = "ExpOrGammaOrChisq"
          )

## Class: Chi squared distribution
setClass("Chisq",
          prototype = prototype(
                      r = function(n){ rchisq(n, df = 1, ncp = 0) },
                      d = function(x, log = FALSE)
                                  { dchisq(x, df = 1, ncp = 0, log = log) },
                      p = function(q, lower.tail = TRUE, log.p = FALSE )
                                  { pchisq(q, df = 1, ncp = 0,
                                           lower.tail = lower.tail,
                                           log.p = log.p) },
                      q = function(p, lower.tail = TRUE, log.p = FALSE )
                                  { qchisq(p, df = 1, ncp = 0,
                                           lower.tail = lower.tail,
                                           log.p = log.p) },
                      param = new("ChisqParameter"),
                     .logExact = TRUE,
                     .lowerExact = TRUE
                      ),
          contains = "ExpOrGammaOrChisq"
          )

## Class: Laplace or Double Exponential distribution
setClass("DExp",
          prototype = prototype(
                      r = function(n){
                              (2*rbinom(n ,size = 1, prob = 0.5)-1) *
                                 rexp(n, rate = 1)
                              },
                      d = function(x, log = FALSE)
                              { d0 <-  dexp(abs(x), rate = 1, log = log)
                                d0 <- if (log) d0-log(2) else d0 <- d0 / 2
                                return(d0) },
                      p = function(q, lower.tail = TRUE, log.p = FALSE ){
                                   if (!lower.tail) q <- -q
                                   p0 <- ifelse( q <= 0,
                                                 0.5*pexp(-q, rate = 1,
                                                          lower.tail = FALSE),
                                                 0.5 + 0.5*pexp( q, rate = 1)
                                                    )
                                   if (log.p)  p0 <- log(p0)
                                   return(p0)
                                   },
                      q = function(p, lower.tail = TRUE, log.p = FALSE ){
                              if (log.p) p <- exp(p)
                              if (!lower.tail) p <- 1-p
                              ifelse( p <= 0.25,
                                  -qexp(2*p, rate = 1, lower.tail = FALSE),
                                  ifelse( p <= 0.5,
                                      -qexp(1-2*p, rate = 1),
                                      ifelse( p <= 0.75   ,
                                          qexp(2*p - 1, rate = 1),
                                          qexp(2*(1-p), rate = 1,
                                               lower.tail = FALSE)
                                            )
                                         )
                                     )},
                      param = new("ExpParameter", name =
                      gettext("Parameter of a Laplace/Double Exponential distribution")
                                 ),
                     .logExact = TRUE,
                     .lowerExact = TRUE,
                     Symmetry = new("SphericalSymmetry", 
                                     type = "univariate symmetric distribution",
                                     SymmCenter = 0)
                      ),
          contains = "AbscontDistribution"
          )

## Class: CauchyDistribution
setClass("Cauchy",
          prototype = prototype(
                      r = function(n){ rcauchy(n,  location = 0, scale = 1) },
                      d = function(x, log = FALSE){
                              dcauchy(x,  location = 0, scale = 1, log = log)
                                          },
                      p = function(q, lower.tail = TRUE, log.p = FALSE ){
                              pcauchy(q,  location = 0, scale = 1,
                                      lower.tail = lower.tail, log.p = log.p)
                                          },
                      q = function(p, lower.tail = TRUE, log.p = FALSE ){
                              qcauchy(p,  location = 0, scale = 1,
                                      lower.tail = lower.tail, log.p = log.p)
                                          },
                      param = new("CauchyParameter"),
                     .logExact = TRUE,
                     .lowerExact = TRUE,
                     Symmetry = new("SphericalSymmetry", 
                                     type = "univariate symmetric distribution",
                                     SymmCenter = 0)
                      ),
          contains = "AbscontDistribution"
          )

## Class: normal distribution
setClass("Norm",
          prototype = prototype(
                      r = function(n){ rnorm(n, mean = 0, sd = 1) },
                      d = function(x, log = FALSE)
                              { dnorm(x, mean = 0, sd = 1, log = log) },
                      p = function(q, lower.tail = TRUE, log.p = FALSE )
                              { pnorm(q, mean = 0, sd = 1,
                                      lower.tail = lower.tail, log.p = log.p) },
                      q = function(p, lower.tail = TRUE, log.p = FALSE )
                              { qnorm(p, mean = 0, sd = 1,
                                      lower.tail = lower.tail, log.p = log.p) },
                      param = new("UniNormParameter"),
                     .logExact = TRUE,
                     .lowerExact = TRUE,
                     Symmetry = new("SphericalSymmetry", 
                                     type = "univariate symmetric distribution",
                                     SymmCenter = 0)
                      ),
          contains = "AbscontDistribution"
          )

## Class: lognormal distribution
setClass("Lnorm",
          prototype = prototype(
                      r = function(n){ rlnorm(n, meanlog = 0, sdlog = 1) },
                      d = function(x, log = FALSE){
                              dlnorm(x, meanlog = 0, sdlog = 1, log = log)
                                          },
                      p = function(q, lower.tail = TRUE, log.p = FALSE ){
                              plnorm(q, meanlog = 0, sdlog = 1,
                                     lower.tail = lower.tail, log.p = log.p)
                                          },
                      q = function(p, lower.tail = TRUE, log.p = FALSE ){
                              qlnorm(p, meanlog = 0, sdlog = 1,
                                     lower.tail = lower.tail, log.p = log.p)
                                          },
                      param = new("LnormParameter"),
                     .logExact = TRUE,
                     .lowerExact = TRUE
                      ),
          contains = "AbscontDistribution"
          )

## Class: uniform distribution
setClass("Unif",
          prototype = prototype(
                      r = function(n){ runif(n, min = 0, max = 1) },
                      d = function(x, log = FALSE)
                              { dunif(x,  min = 0, max = 1, log = log) },
                      p = function(q, lower.tail = TRUE, log.p = FALSE )
                              { punif(q,  min = 0, max = 1,
                                      lower.tail = lower.tail, log.p = log.p) },
                      q = function(p, lower.tail = TRUE, log.p = FALSE )
                              { qunif(p,  min = 0, max = 1,
                                      lower.tail = lower.tail, log.p = log.p) },
                      param = new("UnifParameter"),
                     .logExact = TRUE,
                     .lowerExact = TRUE,
                     Symmetry = new("SphericalSymmetry", 
                                     type = "univariate symmetric distribution",
                                     SymmCenter = .5)
                      ),
          contains = "AbscontDistribution"
          )

## Class: F distribution
setClass("Fd",
          prototype = prototype(
                      r = function(n){ rf(n,  df1 = 1, df2 = 1, ncp = 0) },
                      d = function(x, log = FALSE){
                              df(x, df1 = 1, df2 = 1, ncp = 0, log = log)
                                          },
                      p = function(q, lower.tail = TRUE, log.p = FALSE ){
                              pf(q, df1 = 1, df2 = 1, ncp = 0,
                                 lower.tail = lower.tail, log.p = log.p)
                                          },
                      q = function(p, lower.tail = TRUE, log.p = FALSE ){
                              qf(p, df1 = 1, df2 = 1, ncp = 0,
                                 lower.tail = lower.tail, log.p = log.p)
                                          },
                      param = new("FParameter"),
                     .logExact = TRUE,
                     .lowerExact = TRUE
                      ),
          contains = "AbscontDistribution"
          )

## Class: Student distribution
setClass("Td",
          prototype = prototype(
                      r = function(n){ rt(n,  df = 1, ncp = 0) },
                      d = function(x, log = FALSE)
                              { dt(x,  df = 1, ncp = 0, log = log) },
                      p = function(q, lower.tail = TRUE, log.p = FALSE )
                              { pt(q,  df = 1, ncp = 0,
                                   lower.tail = lower.tail, log.p = log.p) },
                      q = function(p, lower.tail = TRUE, log.p = FALSE )
                              { qt(p,  df = 1, ncp = 0,
                                   lower.tail = lower.tail, log.p = log.p) },
                      param = new("TParameter"),
                     .logExact = TRUE,
                     .lowerExact = TRUE,
                     Symmetry = new("SphericalSymmetry", 
                                     type = "univariate symmetric distribution",
                                     SymmCenter = 0)
                      ),
          contains = "AbscontDistribution"
          )


## Class: logistic distribution
setClass("Logis",
          prototype = prototype(
                      r = function(n){ rlogis(n, location = 0, scale = 1) },
                      d = function(x, log = FALSE){
                              dlogis(x, location = 0, scale = 1, log = log)
                                          },
                      p = function(q, lower.tail = TRUE, log.p = FALSE ){
                              plogis(q, location = 0, scale = 1,
                                     lower.tail = lower.tail, log.p = log.p)
                                          },
                      q = function(p, lower.tail = TRUE, log.p = FALSE ){
                              qlogis(p, location = 0, scale = 1,
                                     lower.tail = lower.tail, log.p = log.p)
                                          },
                      param = new("LogisParameter"),
                     .logExact = TRUE,
                     .lowerExact = TRUE
                      ),
          contains = "AbscontDistribution"
          )

## Class: BetaDistribution
setClass("Beta",
          prototype = prototype(
                      r = function(n){
                              rbeta(n,  shape1 = 1, shape2 = 1, ncp = 0)
                                     },
                      d = function(x, log = FALSE){
                              dbeta(x,  shape1 = 1, shape2 = 1, ncp = 0,
                                    log = log)
                                          },
                      p = function(q, lower.tail = TRUE, log.p = FALSE ){
                              pbeta(q,  shape1 = 1, shape2 = 1, ncp = 0,
                                    lower.tail = lower.tail, log.p = log.p)
                                          },
                      q = function(p, lower.tail = TRUE, log.p = FALSE ){
                              qbeta(p,  shape1 = 1, shape2 = 1, ncp = 0,
                                    lower.tail = lower.tail, log.p = log.p)
                                          },
                      param = new("BetaParameter"),
                     .logExact = TRUE,
                     .lowerExact = TRUE
                      ),
          contains = "AbscontDistribution"
          )

## Class: Weibull distribution
setClass("Weibull",
          prototype = prototype(
                      r = function(n){ rweibull(n, shape = 1, scale = 1) },
                      d = function(x, log = FALSE){
                              dweibull(x, shape = 1, scale = 1, log = log)
                                          },
                      p = function(q, lower.tail = TRUE, log.p = FALSE ){
                              pweibull(q, shape = 1, scale = 1,
                                       lower.tail = lower.tail, log.p = log.p)
                                          },
                      q = function(p, lower.tail = TRUE, log.p = FALSE ){
                              qweibull(p, shape = 1, scale = 1,
                                       lower.tail = lower.tail, log.p = log.p)
                                          },
                      param = new("WeibullParameter"),
                     .logExact = TRUE,
                     .lowerExact = TRUE
                      ),
          contains = "AbscontDistribution"
          )

## Class: Arcsine distribution
setClass("Arcsine",
          prototype = prototype(
                      r = function(n){ sin((runif(n)-.5)*pi) },
                      d = function(x, log = FALSE){
                              x0 <- (abs(x)<1-.Machine$double.eps)
                              x1 <- x^2*x0
                              d <-  x0/sqrt(1-x1)/pi
                              d[.isEqual(abs(x),1)] <- Inf
                              if(log) d<- log(d)
                              return(d)},
                      p = function(q, lower.tail = TRUE, log.p = FALSE ){
                              if(!lower.tail) q<- -q
                              q <- pmin(pmax(q,-1),1)
                              p <- asin(q)/pi+1/2
                              if(log.p) p <- log(p)
                              return(p)},
                      q = function(p, lower.tail = TRUE, log.p = FALSE ){
                              if(log.p) p <- exp(p)
                              p1 <- p
                              p1[p<0|p>1] <- 0.5
                              if(!lower.tail) p1 <- 1-p1
                              q <- sin( (p1-1/2)*pi)
                              q[p<0|p>1] <- NA
                              q[.isEqual(p,0)] <- -1
                              q[.isEqual(p,1)] <-  1
                              return(q)},
                     .logExact = TRUE,
                     .lowerExact = TRUE,
                     Symmetry = new("SphericalSymmetry", 
                                     type = "univariate symmetric distribution",
                                     SymmCenter = 0)
                      ),
          contains = "AbscontDistribution"
          )


## inbetween-Class: AffLinAbscontDistribution

setClass("AffLinAbscontDistribution", 
          representation = representation(a = "numeric", b = "numeric",
          X0 = "AbscontDistribution"),
          prototype = prototype(a = 1, b = 0, X0 = new("Norm")),
          contains = "AbscontDistribution"
          )

### ---- discrete distributions ---- ###

## DiscreteDistribution
setClass("DiscreteDistribution", 
          representation = representation(support = "numeric"),
          prototype = prototype(
                      r = function(n){ rbinom(n, size=1, prob=0.5) },
                      d = function(x, log = FALSE)
                              { dbinom(x, size=1, prob=0.5, log = log) },
                      p = function(q, lower.tail = TRUE, log.p = FALSE )
                              { pbinom(q, size=1, prob=0.5, 
                                       lower.tail = lower.tail, log.p = log.p) },
                      q = function(p, lower.tail = TRUE, log.p = FALSE )
                              { qbinom(p, size=1, prob=0.5, 
                                       lower.tail = lower.tail, log.p = log.p) },
                      img = new("Reals"),
                      support = 0:1 
                      ), 
          contains = "UnivariateDistribution"
          )

## LatticeDistribution
setClass("LatticeDistribution", 
          representation = representation(lattice = "Lattice"),
          prototype = prototype(lattice = new("Lattice")),
          contains = "DiscreteDistribution"
          )


### Class: Dirac distribution
setClass("Dirac",
          prototype = prototype(
                      r = function(n){ array(0, n)},
                      d = function(x, log)
                              { y <- rep(0,length(x))
                                d0 <- as.numeric(x == y)
                                if(log) d0 <- log(d0)
                                return(d0)
                              },
                      p = function(q, lower.tail = TRUE, log.p = FALSE )
                              { p0 <- as.numeric(q + 10^-10 >= 0)
                                if (!lower.tail) p0 <- 1-p0
                                if (log.p) p0 <- log(p0)
                                return(p0)
                              },
                      q = function(p, lower.tail = TRUE, log.p = FALSE )
                             {  if (log.p) p <- exp(p)
                                if(any((p < 0)|(p > 1)))
                                   warning("q Method of class Dirac produced NaN's.")
                                q0 <- 0 * p
                                q0[(p<0) | (p>1)] <- NaN
                                return(q0)
                              },
                      param = new("DiracParameter"),
                      support = 0,
                      lattice = new("Lattice",
                                pivot = 0, width = 1, Length = 1, name =
                                gettext("lattice of a Dirac distribution")
                                ),
                     .logExact = TRUE,
                     .lowerExact = TRUE
                      ),
          contains = "LatticeDistribution"
          )

## Class: Poisson distribution
setClass("Pois",
          prototype = prototype(
                      r = function(n){ rpois(n, lambda = 1) },
                      d = function(x, log = FALSE)
                              { dpois(x, lambda = 1, log = log) },
                      p = function(q, lower.tail = TRUE, log.p = FALSE )
                              { ppois(q, lambda = 1, lower.tail = lower.tail,
                                      log.p = log.p) },
                      q = function(p, lower.tail = TRUE, log.p = FALSE )
                              { qpois(p, lambda = 1, lower.tail = lower.tail,
                                      log.p = log.p) },
                      img = new("Naturals"),
                      param = new("PoisParameter"),
                      support = seq( 0,
                                     qpois(getdistrOption("TruncQuantile"),
                                           lambda = 1, lower.tail = FALSE),
                                     by = 1
                                    ),
                      lattice = new("Lattice",
                                pivot = 0, width = 1, Length = Inf, name =
                                gettext("lattice of a Poisson distribution")
                                ),
                     .logExact = TRUE,
                     .lowerExact = TRUE
                      ),
          contains = "LatticeDistribution"
          )

## Class: negative binomial distribution
setClass("Nbinom",
          prototype = prototype(
                      r = function(n){ rnbinom(n, size = 1, prob = 0.5) },
                      d = function(x, log = FALSE){
                              dnbinom(x, size = 1, prob = 0.5, log = log)
                                          },
                      p = function(q, lower.tail = TRUE, log.p = FALSE ){
                              pnbinom(q, size = 1, prob = 0.5,
                                      lower.tail = lower.tail, log.p = log.p)
                                          },
                      q = function(p, lower.tail = TRUE, log.p = FALSE ){
                              qnbinom(p, size = 1, prob = 0.5,
                                      lower.tail = lower.tail, log.p = log.p)
                                          },
                      img = new("Naturals"),
                      param = new("NbinomParameter"),
                      support = seq( 0,
                                     qnbinom(
                                        getdistrOption("TruncQuantile"),
                                        size = 1, prob = 0.5, lower.tail = FALSE
                                            ),
                                     by = 1
                                     ),
                      lattice = new("Lattice",
                                pivot = 0, width = 1, Length = Inf, name =
                                gettext(
                                  "lattice of a Negative Binomial distribution"
                                       )
                                ),
                     .logExact = TRUE,
                     .lowerExact = TRUE
                      ),
          contains = "LatticeDistribution"
          )

## Class: hypergeometric distribution
setClass("Hyper",
          prototype = prototype(
                      r = function(nn){ rhyper(nn, m = 1, n = 1, k = 1) },
                      d = function(x, log = FALSE){
                              dhyper(x, m = 1, n = 1, k = 1, log = log)
                                          },
                      p = function(q, lower.tail = TRUE, log.p = FALSE ){
                              phyper(q, m = 1, n = 1, k = 1,
                                     lower.tail = lower.tail, log.p = log.p)
                                          },
                      q = function(p, lower.tail = TRUE, log.p = FALSE ){
                              qhyper(p, m = 1, n = 1, k = 1,
                                     lower.tail = lower.tail, log.p = log.p)
                                          },
                      img = new("Naturals"),
                      param = new("HyperParameter"),
                      support = 0:1,
                      lattice = new("Lattice",
                                pivot = 0, width = 1, Length = 2, name =
                                gettext(
                                  "lattice of a Hypergeometric distribution"
                                       )
                                ),
                     .logExact = TRUE,
                     .lowerExact = TRUE
                      ),
          contains = "LatticeDistribution"
          )

## Class: binomial distribution
setClass("Binom",
          prototype = prototype(
                      r = function(n){ rbinom(n, size = 1,prob = 0.5) },
                      d = function(x, log = FALSE){
                              dbinom(x, size = 1, prob = 0.5, log = log)
                                          },
                      p = function(q, lower.tail = TRUE, log.p = FALSE ){
                              pbinom(q, size = 1, prob = 0.5,
                                     lower.tail = lower.tail, log.p = log.p)
                                          },
                      q = function(p, lower.tail = TRUE, log.p = FALSE ){
                              qbinom(p, size = 1, prob = 0.5,
                                     lower.tail = lower.tail, log.p = log.p)
                                          },
                      img = new("Naturals"),
                      param = new("BinomParameter"),
                      support = 0:1,
                      lattice = new("Lattice",
                                pivot = 0, width = 1, Length = 2, name =
                                gettext(
                                  "lattice of a Binomial distribution"
                                       )
                                ),
                     .logExact = TRUE,
                     .lowerExact = TRUE
                      ),
          contains = "LatticeDistribution"
          )

## Class: geometric distribution
setClass("Geom",
          prototype = prototype(
                      r = function(n){ rgeom(n, prob = 0.5) },
                      d = function(x, log = FALSE)
                              { dgeom(x, prob = 0.5, log = log) },
                      p = function(q, lower.tail = TRUE, log.p = FALSE )
                              { pgeom(q, prob = 0.5, lower.tail = lower.tail,
                                      log.p = log.p) },
                      q = function(p, lower.tail = TRUE, log.p = FALSE )
                              { qgeom(p, prob = 0.5, lower.tail = lower.tail,
                                      log.p = log.p) },
                      param = new("NbinomParameter", name =
                              gettext("Parameter of a Geometric distribution")
                              ),
                      support = seq( 0,
                                     qgeom(getdistrOption("TruncQuantile"),
                                           prob = 0.5, lower.tail = FALSE),
                                     by = 1
                                    ),
                      lattice = new("Lattice",
                                pivot = 0, width = 1, Length = Inf, name =
                                gettext(
                                  "lattice of a Geometric distribution"
                                       )
                                ),
                     .logExact = TRUE,
                     .lowerExact = TRUE
                      ),
          contains = "Nbinom"
          )

### ---- List of distributions by M. Kohl ---- ###

setClass(Class = "DistrList", 
            prototype = prototype(list(new("Norm"))), 
            contains = "list",
            validity = function(object){
                nrvalues <- length(object)
                for(i in 1:nrvalues)
                    if(!is(object[[i]], "Distribution")) 
                        stop("element ", i, " is no 'Distribution'")
                return(TRUE) 
            })

## inbetween-Classes: AffLinDiscreteDistribution, AffLinLatticeDistribution

setClass("AffLinDiscreteDistribution", 
          representation = representation(a = "numeric", b = "numeric",
          X0 = "DiscreteDistribution"),
          prototype = prototype(a = 1, b = 0, X0 = new("Binom")),
          contains = "DiscreteDistribution"
          )

setClass("AffLinLatticeDistribution", 
          contains = c("LatticeDistribution", "AffLinDiscreteDistribution")
          )




################################
##
## Distribution List classes 
##
################################

setClass("UnivarDistrList", 
            prototype = prototype(list(new("Norm"))), 
            contains = "DistrList", 
            validity = function(object){
                nrvalues <- length(object)
                for(i in 1:nrvalues)
                    if(!is(object[[i]], "UnivariateDistribution"))
                        stop("element ", i, " is no 'UnivariateDistribution'")
                return(TRUE) 
            })


#### new from version 2.0: Mixing Distributions

################################
##
## Mixing Distribution classes 
##
################################

setClass("UnivarMixingDistribution",
            representation = representation(mixCoeff = "numeric",
                             mixDistr = "UnivarDistrList",
                             gaps = "OptionalMatrix",
                             support = "numeric",
                             Symmetry = "DistributionSymmetry",
                             .withArith = "logical",
                             .withSim = "logical",
                             .logExact = "logical",
                             .lowerExact = "logical"
                             ),
            prototype = prototype(mixCoeff = 1, 
                                  mixDistr = new("UnivarDistrList"),
                                  gaps = NULL,
                                  support = numeric(0),
                                  Symmetry = new("NoSymmetry"),
                                 .withArith = FALSE,
                                 .withSim = FALSE,
                                 .logExact = TRUE,
                                 .lowerExact = TRUE
                                  ),
            contains = "UnivariateDistribution",
            validity = function(object){
                if(any(object@mixCoeff< -.Machine$double.eps) || 
                   sum(object@mixCoeff)>1+.Machine$double.eps)
                   stop("mixing coefficients are no probabilities")
                return(TRUE)
            })

setClass("UnivarLebDecDistribution",
            representation = representation(mixCoeff = "numeric",
                             mixDistr = "UnivarDistrList"),
            prototype = prototype(mixCoeff = c("acWeight" = 1, 
                                               "discreteWeight" = 0),
                                  mixDistr = new("UnivarDistrList",
                                              list("acPart" = new("Norm"),
                                                   "discretePart" = new("Dirac")
                                                   )
                                  )),
            contains = "UnivarMixingDistribution",
            validity = function(object){
                if (length(object@mixCoeff)!=2)
                    stop("number of mixing components is not 2")
                if (!is(object@mixDistr[[1]], "AbscontDistribution"))
                    stop("first component must be absolutely continuous")
                if (!is(object@mixDistr[[2]], "DiscreteDistribution"))
                    stop("second component must be discrete")
                return(TRUE)
            })

setClass("AffLinUnivarLebDecDistribution",
          representation = representation(a = "numeric", b = "numeric",
          X0 = "UnivarLebDecDistribution"),
          prototype = prototype(a = 1, b = 0, 
                                X0 = new("UnivarLebDecDistribution")),
          contains = "UnivarLebDecDistribution"
          )

         
setClassUnion("UnivDistrListOrDistribution",
               c("UnivarDistrList","UnivariateDistribution"))

setClass("CompoundDistribution", representation=representation(
             NumbOfSummandsDistr = "DiscreteDistribution",
             SummandsDistr = "UnivDistrListOrDistribution"),
          prototype=prototype(NumbOfSummandsDistr = new("Pois"),
              SummandsDistr=new("Norm")),
          contains = "UnivarMixingDistribution"
         )

################################
##
## virtual Distribution class Unions 
##
################################

setClassUnion("AcDcLcDistribution", c("AbscontDistribution",
               "DiscreteDistribution", "UnivarLebDecDistribution",
               "CompoundDistribution"))

setClassUnion("AffLinDistribution", c("AffLinAbscontDistribution",
               "AffLinDiscreteDistribution", "AffLinUnivarLebDecDistribution"))

