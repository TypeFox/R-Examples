.onLoad <- function(lib, pkg){
}

# parameter of Gumbel distribution
setClass("GumbelParameter", representation(loc = "numeric",
                                           scale = "numeric"),
            prototype(name = gettext("parameter of a Gumbel distribution"),
                      loc = 0, scale = 1),
            contains = "Parameter",
            validity = function(object){
                if(length(object@scale) != 1)
                    stop("length of 'scale' is not equal to 1")
                if(length(object@loc) != 1)
                    stop("length of 'loc' is not equal to 1")
                if(object@scale <= 0)
                    stop("'scale' has to be positive")
                else return(TRUE)
            })

# Gumbel distribution
setClass("Gumbel",
            prototype = prototype(r = function(n){ rgumbel(n, loc = 0, scale = 1) },
                                  d = function(x, log){ dgumbel(x, loc = 0, scale = 1, log = FALSE) },
                                  p = function(q, lower.tail = TRUE, log.p = FALSE){
                                         p0 <- pgumbel(q, loc = 0, scale = 1, lower.tail = lower.tail)
                                         if(log.p) return(log(p0)) else return(p0)
                                  },
                                  q = function(p, loc = 0, scale = 1, lower.tail = TRUE, log.p = FALSE){
                                      ## P.R.: changed to vectorized form
                                      p1 <- if(log.p) exp(p) else p

                                      in01 <- (p1>1 | p1<0)
                                      i01 <- distr:::.isEqual01(p1)
                                      i0 <- (i01 & p1<1)
                                      i1 <- (i01 & p1>0)
                                      ii01 <- distr:::.isEqual01(p1) | in01

                                      p0 <- p
                                      p0[ii01] <- if(log.p) log(0.5) else 0.5

                                      q1 <- qgumbel(p0, loc = 0, scale = 1,
                                                    lower.tail = lower.tail)
                                      q1[i0] <- if(lower.tail) -Inf else Inf
                                      q1[i1] <- if(!lower.tail) -Inf else Inf
                                      q1[in01] <- NaN

                                      return(q1)
                                      },
                                  img = new("Reals"),
                                  param = new("GumbelParameter"),
                                  .logExact = FALSE,
                                  .lowerExact = TRUE),
            contains = "AbscontDistribution")

# symmetry of functions
setClass("FunctionSymmetry", contains = c("Symmetry", "VIRTUAL"))

# non-symmetric functions
setClass("NonSymmetric", contains = "FunctionSymmetry", 
            prototype = prototype(type = "non-symmetric function",
                                  SymmCenter = NULL))

# even functions
setClass("EvenSymmetric", contains = "FunctionSymmetry", 
            prototype = prototype(type = "even function",
                                  SymmCenter = numeric(0)))

# odd functions
setClass("OddSymmetric", contains = "FunctionSymmetry", 
            prototype = prototype(type = "odd function",
                                  SymmCenter = numeric(0)))

# list of symmetry types
setClass(Class = "FunSymmList", 
            prototype = prototype(list(new("NonSymmetric"))), 
            contains = "list",
            validity = function(object){
                nrvalues <- length(object)
                for(i in 1:nrvalues)
                    if(!is(object[[i]], "FunctionSymmetry")) 
                        stop("element ", i, " is no 'FunctionSymmetry'")
                return(TRUE) 
            })

# Parameter of a parametric family of probability measures
setClass("ParamFamParameter", 
            representation(main = "numeric",
                           nuisance = "OptionalNumeric",
                           trafo = "matrix"),
            prototype(name = "parameter of a parametric family of probability measures",
                      main = numeric(0), nuisance = NULL, trafo = new("matrix")),
            contains = "Parameter",
            validity = function(object){
                dimension <- length(object@main) + length(object@nuisance)
                if(ncol(object@trafo) != dimension)
                    stop("invalid transformation:\n", 
                         "number of columns of 'trafo' not equal to ", 
                         "dimension of the parameter")
                if(nrow(object@trafo) > dimension)
                    stop("invalid transformation:\n",
                         "number of rows of 'trafo' larger than ", 
                         "dimension of the parameter")
                if(any(!is.finite(object@trafo)))
                    stop("infinite or missing values in 'trafo'")
                return(TRUE)
            })

# family of probability measures
setClass("ProbFamily", representation(name = "character", 
                                      distribution = "Distribution",
                                      distrSymm = "DistributionSymmetry",
                                      props = "character"), 
                       contains = "VIRTUAL")

# parametric family of probability measures
setClass("ParamFamily", representation(param = "ParamFamParameter"), 
            prototype(name = "parametric family of probability measures",
                      distribution = new("Norm"),
                      distrSymm = new("NoSymmetry"),
                      props = character(0),
                      param = new("ParamFamParameter", main = 0, trafo = as.matrix(1))),
            contains = "ProbFamily")

# L_2 differentiable parametric family
setClass("L2ParamFamily", 
            representation(L2deriv = "EuclRandVarList",
                           L2derivSymm = "FunSymmList",
                           L2derivDistr = "DistrList",
                           L2derivDistrSymm = "DistrSymmList",
                           FisherInfo = "PosDefSymmMatrix"), 
            prototype(name = "L_2 differentiable parametric family of probability measures",
                      distribution = new("Norm"),
                      distrSymm = new("NoSymmetry"),
                      param = new("ParamFamParameter", main = 0, trafo = matrix(1)),
                      props = character(0),
                      L2deriv = EuclRandVarList(RealRandVariable(Map = list(function(x){x}), 
                                                                 Domain = Reals())),
                      L2derivSymm = new("FunSymmList"),
                      L2derivDistr = UnivarDistrList(new("Norm")),
                      L2derivDistrSymm = new("DistrSymmList"),
                      FisherInfo = new("PosDefSymmMatrix", matrix(1))),
            contains = "ParamFamily", 
            validity = function(object){
                if(is(object@distribution, "UnivariateCondDistribution"))
                    stop("conditional distributions are not allowed in slot 'distribution'")

                if(!is(object@distrSymm, "NoSymmetry")){
                    if(!is(object@distrSymm@SymmCenter, "numeric"))
                        stop("slot 'SymmCenter' of 'distrSymm' has to be of class 'numeric'")
                    if(length(object@distrSymm@SymmCenter) != dimension(img(object@distribution)))
                        stop("slot 'SymmCenter' of 'distrSymm' has wrong dimension")
                }

                dims <- length(object@param)
                if(ncol(object@FisherInfo) != dims)
                    stop(paste("dimension of 'FisherInfo' should be", dims))
                nrvalues <- numberOfMaps(object@L2deriv)
                if(nrvalues != length(object@L2derivSymm))
                    stop("number of Maps of 'L2deriv' != length of 'L2derivSymm'")
                if(nrvalues != length(object@L2derivDistr))
                    stop("number of Maps of 'L2deriv' != length of 'L2derivDistr'")
                if(nrvalues != length(object@L2derivDistrSymm))
                    stop("number of Maps of 'L2deriv' != length of 'L2derivDistrSymm'")
                if(dimension(Domain(object@L2deriv[[1]])) != dimension(img(object@distribution)))
                    stop("dimension of 'Domain' of 'L2deriv' != dimension of 'img' of 'distribution'")
                if(dimension(object@L2deriv) != dims)
                    stop("dimension of 'L2deriv' != dimension of parameters")
                    
                return(TRUE) 
            })

# neighborhood
setClass("Neighborhood", 
            representation(type = "character",
                           radius = "numeric"), 
            contains = "VIRTUAL")
# unconditional (errors-in-variables) neighborhood
setClass("UncondNeighborhood", contains = c("Neighborhood", "VIRTUAL"))
# unconditional convex contamination neighborhood
setClass("ContNeighborhood", contains = "UncondNeighborhood", 
            prototype = prototype(type = "(uncond.) convex contamination neighborhood",
                                  radius = 0))
# unconditional total variation neighborhood
setClass("TotalVarNeighborhood", contains = "UncondNeighborhood",
            prototype = prototype(type = "(uncond.) total variation neighborhood",
                                  radius = 0))
# robust model
setClass("RobModel", 
            representation(center = "ProbFamily",
                           neighbor = "Neighborhood"), 
            contains = "VIRTUAL")
# robust model with fixed (unconditional) neighborhood
setClass("FixRobModel", 
            prototype = prototype(center = new("ParamFamily"),
                                  neighbor = new("ContNeighborhood")),
            contains = "RobModel",
            validity = function(object){
                if(!is(object@neighbor, "UncondNeighborhood"))
                    stop("'neighbor' is no unconditional neighborhood")
                if(any(object@neighbor@radius < 0 || object@neighbor@radius > 1))
                    stop("neighborhood radius has to be in [0, 1]")
                else return(TRUE)
            })
# robust model with infinitesimal (unconditional) neighborhood
setClass("InfRobModel",  
            prototype = prototype(center = new("L2ParamFamily"),
                                  neighbor = new("ContNeighborhood")),
            contains = "RobModel", 
            validity = function(object){
                if(!is(object@neighbor, "UncondNeighborhood"))
                    stop("'neighbor' is no unconditional neighborhood")
                if(any(object@neighbor@radius < 0))
                    stop("'radius' has to be in [0, Inf]")
                else return(TRUE)
            })
# risks (e.g., risk of estimator)
setClass("RiskType", representation(type = "character"), contains = "VIRTUAL")
# asymptotic risk
setClass("asRisk", contains = c("RiskType", "VIRTUAL"))
# asymptotic covariance
setClass("asCov", contains = "asRisk", 
            prototype = prototype(type = "asymptotic covariance"))
# trace of asymptotic covariance
setClass("trAsCov", contains = "asRisk", 
            prototype = prototype(type = "trace of asymptotic covariance"))
# asymptotic Hampel risk
setClass("asHampel", representation(bound = "numeric"), 
            prototype = prototype(bound = Inf, 
                             type = "trace of asymptotic covariance for given bias bound"),
            contains = "asRisk", 
            validity = function(object){
                if(any(object@bound <= 0))
                    stop("'bound' has to be positive")
                else TRUE
            })
# asymptotic bias
setClass("asBias", contains = "asRisk", 
            prototype = prototype(type = "asymptotic bias"))

# convex asymptotic risk
setClass("asGRisk", contains = c("asRisk", "VIRTUAL")) 
# asymptotic mean square error
setClass("asMSE", contains = "asGRisk", 
            prototype = prototype(type = "asymptotic mean square error"))
# asymptotic under-/overshoot probability
setClass("asUnOvShoot", representation(width = "numeric"), 
            prototype = prototype(type = "asymptotic under-/overshoot probability"),
            contains = "asGRisk",
            validity = function(object){
                if(length(object@width) != 1)
                    stop("length of 'width' has to be 1")
                if(any(object@width <= 0))
                    stop("'width' has to be positive")
                else TRUE
            })
# finite-sample risk
setClass("fiRisk", contains = c("RiskType", "VIRTUAL"))
# finite-sample covariance
setClass("fiCov", contains = "fiRisk", 
            prototype = prototype(type = "finite-sample covariance"))
# trace of finite-sample covariance
setClass("trFiCov", contains = "fiRisk", 
            prototype = prototype(type = "trace of finite-sample covariance"))
# finite-sample Hampel risk
setClass("fiHampel", representation(bound = "numeric"),
            prototype = prototype(bound = Inf, 
                             type = "finite-sample variance for given bias bound"),
            contains = "fiRisk", 
            validity = function(object){
                if(any(object@bound <= 0))
                    stop("'bound' has to be positive")
                else TRUE
            })
# finite-sample mean square error
setClass("fiMSE", contains = "fiRisk", 
            prototype = prototype(type = "finite-sample mean square error"))
# finite-sample bias
setClass("fiBias", contains = "fiRisk", 
            prototype = prototype(type = "finite-sample bias"))
# finite-sample under-/overshoot probability
setClass("fiUnOvShoot", representation(width = "numeric"), 
            prototype = prototype(type = "finite-sample under-/overshoot probability"),
            contains = "fiRisk",
            validity = function(object){
                if(length(object@width) != 1)
                    stop("length of 'width' has to be 1")
                if(any(object@width <= 0))
                    stop("'width' has to be positive")
                else TRUE
            })
# Influence curve/function with domain: EuclideanSpace
setClass("InfluenceCurve", 
            representation(name = "character", 
                           Curve = "EuclRandVarList", 
                           Risks = "list",
                           Infos = "matrix"),
            validity = function(object){
                if(!is(Domain(object@Curve[[1]]), "EuclideanSpace"))
                    stop("The domain of 'Curve' has to be a Euclidean space")
                if(!is.character(object@Infos))
                    stop("'Infos' contains no matrix of characters")
                for(char in names(object@Risks))
                    if(!extends(char, "RiskType"))
                        stop(paste(char, "is no valid 'RiskType'"))
                if(ncol(object@Infos)!=2)
                    stop("'Infos' must have two columns")
                else TRUE
            })
# partial incluence curve
setClass("IC", representation(CallL2Fam = "call"),
            prototype(name = "square integrable (partial) influence curve",
                      Curve = EuclRandVarList(RealRandVariable(Map = list(function(x){x}), 
                                               Domain = Reals())),
                      Risks = list(),
                      Infos = matrix(c(character(0),character(0)), ncol=2,
                                dimnames=list(character(0), c("method", "message"))),
                      CallL2Fam = call("L2ParamFamily")),
            contains = "InfluenceCurve",
            validity = function(object){
                L2Fam <- eval(object@CallL2Fam)
                trafo <- L2Fam@param@trafo
                if(nrow(trafo) != dimension(object@Curve))
                    stop("wrong dimension of 'Curve'")
                if(dimension(Domain(L2Fam@L2deriv[[1]])) != dimension(Domain(object@Curve[[1]])))
                    stop("dimension of 'Domain' of 'L2deriv' != dimension of 'Domain' of 'Curve'")

                return(TRUE)
            })
# (partial) influence curve of contamination type
setClass("ContIC", 
            representation(clip = "numeric",
                           cent = "numeric",
                           stand = "matrix",
                           lowerCase = "OptionalNumeric",
                           neighborRadius = "numeric"), 
            prototype(name = "IC of contamination type",
                      Curve = EuclRandVarList(RealRandVariable(Map = list(function(x){x}), 
                                                    Domain = Reals())),
                      Risks = list(),
                      Infos = matrix(c(character(0),character(0)), ncol=2,
                                dimnames=list(character(0), c("method", "message"))),
                      CallL2Fam = call("L2ParamFamily"),
                      clip = Inf, cent = 0, stand = as.matrix(1),
                      lowerCase = NULL,
                      neighborRadius = 0),
            contains = "IC",
            validity = function(object){
                if(any(object@neighborRadius < 0)) # radius vector?!
                    stop("'neighborRadius' has to be in [0, Inf]")
                if(length(object@cent) != nrow(object@stand))
                    stop("length of centering constant != nrow of standardizing matrix")
                if((length(object@clip) != 1) && (length(object@clip) != length(object@Curve)))
                    stop("length of clipping bound != 1 and != length of 'Curve'")
                if(!is.null(object@lowerCase))
                    if(length(object@lowerCase) != nrow(object@stand))
                        stop("length of 'lowerCase' != nrow of standardizing matrix")
                L2Fam <- eval(object@CallL2Fam)
                if(!identical(dim(L2Fam@param@trafo), dim(object@stand)))
                    stop(paste("dimension of 'trafo' of 'param'", 
                                "!= dimension of 'stand'"))
                return(TRUE)
            })
# (partial) influence curve of total variation type
setClass("TotalVarIC", 
            representation(clipLo = "numeric",
                           clipUp = "numeric",
                           stand = "matrix",
                           lowerCase = "OptionalNumeric",
                           neighborRadius = "numeric"), 
            prototype(name = "IC of total variation type",
                      Curve = EuclRandVarList(RealRandVariable(Map = list(function(x){x}), 
                                                               Domain = Reals())),
                      Risks = list(),
                      Infos = matrix(c(character(0),character(0)), ncol=2,
                                dimnames=list(character(0), c("method", "message"))),
                      CallL2Fam = call("L2ParamFamily"),
                      clipLo = -Inf, clipUp = Inf, stand = as.matrix(1),
                      lowerCase = NULL,
                      neighborRadius = 0),
            contains = "IC",
            validity = function(object){
                if(any(object@neighborRadius < 0)) # radius vector?!
                    stop("'neighborRadius' has to be in [0, Inf]")
                if((length(object@clipLo) != 1) && (length(object@clipLo) != length(object@Curve)))
                    stop("length of lower clipping bound != 1 and != length of 'Curve'")
                if((length(object@clipLo) != 1) && (length(object@clipLo) != length(object@Curve)))
                    stop("length of upper clipping bound != 1 and != length of 'Curve'")
                L2Fam <- eval(object@CallL2Fam)
                if(!identical(dim(L2Fam@param@trafo), dim(object@stand)))
                    stop(paste("dimension of 'trafo' of 'param'", 
                                "!= dimension of 'stand'"))
                return(TRUE)
            })
