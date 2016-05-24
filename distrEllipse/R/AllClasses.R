################################################################################
#
# Parameter Classes
# 
################################################################################


################################################################################
# parameter of Elliptical distribution
################################################################################

setClass("EllipticalParameter", representation(loc = "numeric",
                                           scale = "matrix"),
            prototype(name = 
            gettext("parameter of an elliptically symmetric distribution"),
                      loc=c(0,0), scale = diag(2)),
            contains = "Parameter",
            validity = function(object){
               dim0 <- length(object@loc)
               if(!nrow(object@scale)==dim0) stop("wrong dimensions")
               else return(TRUE)
            })

################################################################################
# parameter of MV norm distribution
################################################################################
setClass("MVNormParameter",
            contains = "EllipticalParameter")

################################################################################
# parameter of MV norm distribution
################################################################################
setClass("MVtParameter",
            representation(df = "numeric", ncp = "numeric"),
            prototype(name = gettext("parameter of multivariate t distribution"),
                      ncp = 0, df = 1),
            contains = "EllipticalParameter",
            validity = function(object){
               dim0 <- length(object@loc)
               if(!distr:::.isNatural(object@df)) stop("'df' must be an integer")
               if(!length(object@ncp)==1) stop("wrong dimension for ncp")
               if(!nrow(object@scale)==dim0) stop("wrong dimensions")
               else return(TRUE)
            })

################################################################################
#
# Distribution Classes
# 
################################################################################


################################################################################
# spherically symmetric distributions
################################################################################
setClass("SphericalDistribution",
            representation = representation(radDistr="UnivariateDistribution",
            Symmetry = "EllipticalSymmetry"),
            prototype = prototype(r = function(n) matrix(rnorm(2*n),ncol=2),
                                  d = function(x, log = FALSE){
                                      r2 <- sum(x^2)
                                      lg <- -p/2*log(2*pi)-r2/2;
                                      return(if(log) lg else exp(lg))},
                                  radDistr = sqrt(Chisq(df=2)),
                                  Symmetry = SphericalSymmetry(0)),
            contains = "MultivariateDistribution")

################################################################################
# elliptically symmetric distributions
################################################################################

setClass("EllipticalDistribution",
            prototype = prototype(param = new("EllipticalParameter")),
            contains = "SphericalDistribution")

################################################################################
# Multivariate Normal
################################################################################
setClass("MVNormDistribution",
            prototype = prototype(
            r = function(n){rmvnorm(n, mean = c(0,0), sigma = diag(2))},
            d = function(x, log = FALSE){dmvnorm(x, mean = c(0,0), sigma = diag(2), log = log)},
            p = function(lower=-Inf, upper=Inf){
                  pmvnorm(lower=lower, upper=upper, mean = c(0,0), sigma = diag(2))},
            q = function(p, interval = c(-10, 10), tail = c("lower.tail",
                         "upper.tail", "both.tails")){
                  qmvnorm(p = p, interval = interval, tail = tail,
                          mean = c(0,0), sigma = diag(2))},
            param = new("MVNormParameter", 
               name = gettext("parameter of multivariate normal distribution"))),
            contains = "EllipticalDistribution")

################################################################################
# Multivariate T
################################################################################

setClass("MVtDistribution",
            prototype = prototype(
            r = function(n){rmvt(n)},
            d = function(x, log = FALSE){dmvt(x, delta=0, sigma = diag(2), log = log)},
            p = function(lower=-Inf, upper=Inf){
                  pmvt(lower=lower, upper=upper, delta=0, df=1, sigma = diag(2))},
            q = function(p, interval = c(-10, 10), tail = c("lower.tail",
                         "upper.tail", "both.tails")){
                  qmvt(p = p, interval = interval, tail = tail, df = 1, delta = 0,
                       sigma = diag(2))},
            param = new("MVtParameter")),
            contains = "EllipticalDistribution")

################################
##
## Distribution List classes 
##
################################

setClass("MVDistrList",
            prototype = prototype(list(new("MVNormDistribution"))),
            contains = "DistrList", 
            validity = function(object){
                nrvalues <- length(object)
                dim0 <- object[[1]]@img@dimension
            #    if (dim0 == 1) return(getValidity(getClass("UnivarDistrList"))(object))
                for(i in 1:nrvalues){
                    if(!is(object[[i]], "MultivariateDistribution"))
                        stop("Element ", i, " is no 'MultivariateDistribution'")
                    if(!object[[i]]@img@dimension==dim0)
                        stop("Dimension mismatch in element ", i, ".")
                }
                return(TRUE) 
            })

setClassUnion("MultivarDistrList", c("MVDistrList","UnivarDistrList"))

################################
##
## Mixing Distribution classes 
##
################################

setClass("MultivarMixingDistribution",
            representation = representation(mixCoeff = "numeric",
                             mixDistr = "MultivarDistrList",
                             Symmetry = "DistributionSymmetry",
                             .withArith = "logical",
                             .withSim = "logical",
                             .logExact = "logical",
                             .lowerExact = "logical"
                             ),
            prototype = prototype(mixCoeff = 1, 
                                  mixDistr = new("MVDistrList"),
                                  Symmetry = new("NoSymmetry"),
                                 .withArith = FALSE,
                                 .withSim = FALSE,
                                 .logExact = TRUE,
                                 .lowerExact = TRUE
                                  ),
            contains = "MultivariateDistribution",
            validity = function(object){
                if(any(object@mixCoeff< -.Machine$double.eps) || 
                   sum(object@mixCoeff)>1+.Machine$double.eps)
                   stop("mixing coefficients are no probabilities")
                return(TRUE)
            })
