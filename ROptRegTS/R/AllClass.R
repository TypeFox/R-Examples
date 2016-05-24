.onLoad <- function(lib, pkg){
}

# Regression type families
setClass("RegTypeFamily", 
            representation(ErrorDistr = "Distribution", 
                           ErrorSymm = "DistributionSymmetry",
                           RegDistr = "Distribution",
                           RegSymm = "DistributionSymmetry",
                           Regressor = "EuclRandVariable"),
            prototype(name = "regression type family", 
                      distribution = LMCondDistribution(),
                      distrSymm = new("NoSymmetry"),
                      ErrorDistr = Norm(),
                      ErrorSymm = new("NoSymmetry"),
                      RegDistr = Norm(),
                      RegSymm = new("NoSymmetry"),
                      Regressor = RealRandVariable(Map = list(function(x){x}), Domain = Reals()),
                      param = new("ParamFamParameter", main = 0, trafo = matrix(1)),
                      props = character(0)),
            contains = "ParamFamily",
            validity = function(object){
                if(!is(object@distribution, "UnivariateCondDistribution"))
                    stop("distribution has to be of class 'UnivariateCondDistribution'")
                if(length(object@Regressor) != 1)
                    stop("'Regressor' has to be of length 1")
                if(is(object@ErrorDistr, "UnivariateCondDistribution"))
                    stop("'ErrorDistr' has to be an unconditional distribution")
                if(is(object@RegDistr, "UnivariateCondDistribution"))
                    stop("'RegrDistr' has to be an unconditional distribution")
                if(dimension(Domain(object@Regressor)) != dimension(img(object@RegDistr)))
                    stop("dimension of 'Domain' of 'Regressor' has to be identical to",
                         "dimension of 'img' of 'RegDistr'")
            })
# L2 differentiable regression type model
setClass("L2RegTypeFamily", 
            representation(L2deriv = "EuclRandVarList",
                           ErrorL2deriv = "EuclRandVarList",
                           ErrorL2derivSymm = "FunSymmList", 
                           ErrorL2derivDistr = "DistrList",
                           ErrorL2derivDistrSymm = "DistrSymmList", 
                           FisherInfo = "PosDefSymmMatrix"), 
            prototype(name = "L2 differentiable regression type family", 
                      distribution = LMCondDistribution(),
                      distrSymm = new("NoSymmetry"),
                      ErrorDistr = Norm(),
                      ErrorSymm = new("NoSymmetry"),
                      RegDistr = Norm(),
                      RegSymm = new("NoSymmetry"),
                      Regressor = RealRandVariable(Map = list(function(x){x}), Domain = Reals()),
                      param = new("ParamFamParameter", main = 0, trafo = matrix(1)),
                      props = character(0),
                      L2deriv = EuclRandVarList(RealRandVariable(Map = list(function(x){x[1]*x[2]}), 
                                                            Domain = EuclideanSpace(dimension=2))),
                      ErrorL2deriv = EuclRandVarList(RealRandVariable(Map = list(function(x){x}), Domain = Reals())),
                      ErrorL2derivSymm = new("FunSymmList"),
                      ErrorL2derivDistr = UnivarDistrList(Norm()),
                      ErrorL2derivDistrSymm = new("DistrSymmList"),
                      FisherInfo = new("PosDefSymmMatrix", matrix(1))),
            contains = "RegTypeFamily", 
            validity = function(object){
                if(dimension(Domain(object@ErrorL2deriv[[1]])) != dimension(img(object@ErrorDistr)))
                    stop("dimension of 'Domain' of 'ErrorL2deriv' != ",
                         "dimension of 'img' of 'ErrorDistr'")
                if(dimension(Domain(object@L2deriv[[1]])) != (dimension(img(object@ErrorDistr))
                                                       + dimension(img(object@RegDistr))))
                    stop("'Domain' of 'L2deriv' has wrong dimension")

                dims <- length(object@param)
                if(ncol(object@FisherInfo) != dims)
                    stop(paste("dimension of 'FisherInfo' should be", dims))

                nrvalues <- numberOfMaps(object@ErrorL2deriv)
                if(nrvalues != length(object@ErrorL2derivSymm))
                    stop("number of Maps of 'ErrorL2deriv' != length of 'ErrorL2derivSymm'")
                if(nrvalues != length(object@ErrorL2derivDistr))
                    stop("number of Maps of 'ErrorL2deriv' != length of 'ErrorL2derivDistr'")
                if(nrvalues != length(object@ErrorL2derivDistrSymm))
                    stop("number of Maps of 'ErrorL2deriv' != length of 'ErrorL2derivDistrSymm'")
                if(dimension(Domain(object@ErrorL2deriv[[1]])) != dimension(img(object@ErrorDistr)))
                    stop("dimension of 'Domain' of 'L2deriv' != dimension of 'img' of 'ErrorDistr'")
                if(dimension(object@L2deriv) != dims)
                    stop("dimension of 'L2deriv' != dimension of parameters")

                return(TRUE) 
            })
# conditional (error-free-variables) neighborhood
setClass("CondNeighborhood", 
            representation(radiusCurve = "function"), 
            contains = c("Neighborhood", "VIRTUAL"),
            validity = function(object){
                if(length(formals(object@radiusCurve)) != 1)
                    stop("'radiusCurve' has to be a function of one argument")
                if(names(formals(object@radiusCurve)) != "x")
                    stop("'radiusCurve' has to be a function with argument name = 'x'")
            })
# conditional convex contamination neighborhood
setClass("CondContNeighborhood", 
            prototype = prototype(type = "conditional convex contamination neighborhood",
                                  radius = 0,
                                  radiusCurve = function(x){1}),
            contains = "CondNeighborhood")
# conditional total variaton neighborhood
setClass("CondTotalVarNeighborhood", 
            prototype = prototype(type = "conditional total variation neighborhood",
                                  radius = 0, 
                                  radiusCurve = function(x){1}),
            contains = "CondNeighborhood")
# average conditional neighborhood
setClass("AvCondNeighborhood", representation(exponent = "numeric"),
            contains = c("CondNeighborhood", "VIRTUAL"))
# average conditional neighborhood (exponent = 1)
setClass("Av1CondNeighborhood", 
            contains = c("AvCondNeighborhood", "VIRTUAL"),
            validity = function(object){
                if(object@exponent != 1)
                    stop("exponent has to be 1")
            })
# average conditional convex contamination neighborhood (exponent = 1)
setClass("Av1CondContNeighborhood", 
            prototype = prototype(type = "average conditional convex contamination neighborhood",
                                  radius = 0,
                                  radiusCurve = function(x){1},
                                  exponent = 1),
            contains = c("Av1CondNeighborhood"))
# average conditional total variation neighborhood (exponent = 1)
setClass("Av1CondTotalVarNeighborhood", 
            prototype = prototype(type = "average conditional total variation neighborhood",
                                  radius = 0,
                                  radiusCurve = function(x){1},
                                  exponent = 1),
            contains = c("Av1CondNeighborhood"))
# average square conditional neighborhood (exponent = 2)
setClass("Av2CondNeighborhood", 
            contains = c("AvCondNeighborhood", "VIRTUAL"),
            validity = function(object){
                if(object@exponent != 2)
                    stop("exponent has to be 2")
            })
# average square conditional convex contamination neighborhood (exponent = 2)
setClass("Av2CondContNeighborhood", 
            prototype = prototype(type = "average square conditional convex contamination neighborhood",
                                  radius = 0,
                                  radiusCurve = function(x){1},
                                  exponent = 2),
            contains = c("Av2CondNeighborhood"))
# robust regression-type model with fixed 
# (conditional or unconditional) neighborhood
setClass("FixRobRegTypeModel", 
            prototype = prototype(center = new("RegTypeFamily"),
                                  neighbor = new("ContNeighborhood")),
            contains = "RobModel",
            validity = function(object){
                if(!is(object@center, "RegTypeFamily"))
                    stop("center has to be a regression type family")
                if(any(object@neighbor@radius < 0 || object@neighbor@radius > 1))
                    stop("neighborhood radius has to be in [0, 1]")
                if(is(object@neighbor, "CondNeighborhood")){
                    D1 <- object@center@RegDistr
                    radCurve <- object@neighbor@radiusCurve
                    if(is(D1, "UnivariateDistribution")){
                        if(is(D1, "AbscontDistribution")){
                            xlo <- ifelse(is.finite(q(D1)(0)), q(D1)(0), q(D1)(distr::TruncQuantile))
                            xup <- ifelse(is.finite(q(D1)(1)), q(D1)(1), q(D1)(1 - distr::TruncQuantile))
                            x <- seq(from = xlo, to = xup, by = 1e-3)
                        }else{
                            if(is(Regressor, "DiscreteDistribution"))
                                x <- support(D1) 
                            else
                                x <- unique(r(D1)(1e5))
                        }
                        if(length(radCurve(x[1])) != 1)
                            stop("'radiusCurve' has to be a real-valued function")
                    }else{
                        if(is(D1, "DiscreteMVDistribution"))
                            x <- support(D1)
                        else
                            x <- r(D1)(1e5)
                        if(length(radCurve(x[1,])) != 1)
                            stop("'radiusCurve' has to be a real-valued function")
                    }
                    if(min(radCurve(x)) < 0 || max(radCurve(x)) > 1)
                        stop("'radiusCurve' has to be a non-negative function",
                             " with values in [0,1]")
                    if(!is.finite(E(D1, radCurve)))
                        stop("'radiusCurve' has an infinite integral")
                    if(is(object@neighbor, "AvCondNeighborhood")){
                        alpha <- object@neighbor@exponent
                        radCurve1 <- function(x, radCurve, alpha){ radCurve(x)^alpha }
                        constr <- E(D1, radCurve1, radCurve = radCurve, alpha = alpha)^(1/alpha)
                        if(constr > 1+.Machine$double.eps^0.5)
                            stop("'radiusCurve' does not fulfill the norm constraint")
                    }
                }
                return(TRUE)
            })
# robust regression-type model with infinitesimal
# (conditional or unconditional) neighborhood
setClass("InfRobRegTypeModel", 
            prototype = prototype(center = new("L2RegTypeFamily"),
                                  neighbor = new("ContNeighborhood")),
            contains = "RobModel", 
            validity = function(object){
                if(!is(object@center, "L2RegTypeFamily"))
                    stop("'center' is no 'L2RegTypeFamily'")
                if(any(object@neighbor@radius < 0)) # radius vector?!
                    stop("'radius' has to be in [0, Inf]")
                if(is(object@neighbor, "CondNeighborhood")){
                    D1 <- object@center@RegDistr
                    radCurve <- object@neighbor@radiusCurve
                    if(is(D1, "UnivariateDistribution")){
                        if(is(D1, "AbscontDistribution")){
                            xlo <- ifelse(is.finite(q(D1)(0)), q(D1)(0), q(D1)(distr::TruncQuantile))
                            xup <- ifelse(is.finite(q(D1)(1)), q(D1)(1), q(D1)(1 - distr::TruncQuantile))
                            x <- seq(from = xlo, to = xup, by = 1e-3)
                        }else{
                            if(is(Regressor, "DiscreteDistribution"))
                                x <- support(D1) 
                            else
                                x <- unique(r(D1)(1e5))
                        }
                        if(length(radCurve(x[1])) != 1)
                            stop("'radiusCurve' has to be a real-valued function")
                    }else{
                        if(is(D1, "DiscreteMVDistribution"))
                            x <- support(D1)
                        else
                            x <- r(D1)(1e5)
                        if(length(radCurve(x[1,])) != 1)
                            stop("'radiusCurve' has to be a real-valued function")
                    }
                    if(min(radCurve(x)) < 0)
                        stop("'radiusCurve' has to be a non-negative function")
                    if(!is.finite(E(D1, radCurve)))
                        stop("'radiusCurve' has an infinite integral")
                    if(is(object@neighbor, "AvCondNeighborhood")){
                        alpha <- object@neighbor@exponent
                        radCurve1 <- function(x, radCurve, alpha){ radCurve(x)^alpha }
                        constr <- E(D1, radCurve1, radCurve = radCurve, alpha = alpha)^(1/alpha)
                        if(constr > 1+.Machine$double.eps^0.5)
                            stop("'radiusCurve' does not fulfill the norm constraint")
                    }
                }else 
                    return(TRUE)
            })
# square integrable, conditionally centered (partial) IC
setClass("CondIC", 
            prototype = prototype(name = "square integrable, conditionally centered (partial) IC",
                                  Curve = EuclRandVarList(EuclRandVariable(Map = list(function(x){x[1]*x[2]}),
                                                Domain = EuclideanSpace(dimension=2), Range = Reals())),
                                  Risks = list(),
                                  Infos = matrix(c(character(0),character(0)), ncol=2,
                                              dimnames=list(character(0), c("method", "message"))),
                                  CallL2Fam = call("L2RegTypeFamily")),
            contains = "IC",
            validity = function(object){
                L2Fam <- eval(object@CallL2Fam)
                if(!is(L2Fam, "L2RegTypeFamily"))
                    stop("'CallL2Fam' has to generate an object of class 'L2RegTypeFamily'")

                return(TRUE)
            })
# square integrable, conditionally centered (partial) IC 
# of contamination type
setClass("CondContIC", 
            representation(clip = "RealRandVariable",
                           cent = "EuclRandVarList",
                           stand = "matrix",
                           lowerCase = "OptionalNumeric",
                           neighborRadius = "numeric",
                           neighborRadiusCurve = "function"), 
            prototype(name = "conditionally centered IC for average conditional contamination neighborhoods",
                      Curve = EuclRandVarList(RealRandVariable(Map = list(function(x){x[1]*x[2]}), 
                                                    Domain = EuclideanSpace(dimension = 2))),
                      Risks = list(),
                      Infos = matrix(c(character(0),character(0)), ncol=2,
                                dimnames=list(character(0), c("method", "message"))),
                      CallL2Fam = call("L2RegTypeFamily"),
                      clip = RealRandVariable(Map = list(function(x){ Inf }), 
                                              Domain = EuclideanSpace(dimension = 1)), 
                      cent = EuclRandVarList(RealRandVariable(Map = list(function(x){numeric(length(x))}),
                                                    Domain = EuclideanSpace(dimension = 2))),
                      stand = as.matrix(1),
                      lowerCase = NULL,
                      neighborRadius = 0,
                      neighborRadiusCurve = function(x){1}),
            contains = "CondIC",
            validity = function(object){
                if(any(object@neighborRadius < 0)) # radius vector?!
                    stop("'neighborRadius' has to be in [0, Inf]")
                if(length(formals(object@neighborRadiusCurve)) != 1)
                    stop("'neighborRadiusCurve' has to be a function of one argument")
                if(names(formals(object@neighborRadiusCurve)) != "x")
                    stop("'neighborRadiusCurve' has to be a function with argument name = 'x'")
                if(dimension(object@cent) != nrow(object@stand))
                    stop("dimension of centering function != nrow of standardizing matrix")
                if(dimension(object@clip) != 1)
                    stop("dimension of clipping function has to be 1")
                if(!is.null(object@lowerCase))
                    if(length(object@lowerCase) != nrow(object@stand))
                        stop("length of 'lowerCase' != nrow of standardizing matrix")
                L2Fam <- eval(object@CallL2Fam)
                if(!identical(dim(L2Fam@param@trafo), dim(object@stand)))
                    stop("dimension of 'trafo' of 'param' != dimension of 'stand'")
                
                return(TRUE)
            })
# square integrable, conditionally centered (partial) IC 
# of contamination type
setClass("Av1CondContIC", 
            representation(clip = "numeric",
                           cent = "EuclRandVarList",
                           stand = "matrix",
                           lowerCase = "OptionalNumeric",
                           neighborRadius = "numeric"), 
            prototype(name = "conditionally centered IC for average conditional contamination neighborhoods",
                      Curve = EuclRandVarList(RealRandVariable(Map = list(function(x){x[1]*x[2]}), 
                                                    Domain = EuclideanSpace(dimension = 2))),
                      Risks = list(),
                      Infos = matrix(c(character(0),character(0)), ncol=2,
                                dimnames=list(character(0), c("method", "message"))),
                      CallL2Fam = call("L2RegTypeFamily"),
                      clip = Inf, 
                      cent = EuclRandVarList(RealRandVariable(Map = list(function(x){numeric(length(x))}),
                                                    Domain = EuclideanSpace(dimension = 2))),
                      stand = as.matrix(1),
                      lowerCase = NULL,
                      neighborRadius = 0),
            contains = "CondIC",
            validity = function(object){
                if(any(object@neighborRadius < 0)) # radius vector?!
                    stop("'neighborRadius' has to be in [0, Inf]")
                if(dimension(object@cent) != nrow(object@stand))
                    stop("dimension of centering function != nrow of standardizing matrix")
                if(length(object@clip) != 1)
                    stop("length of clipping bound has to be 1")
                if(!is.null(object@lowerCase))
                    if(length(object@lowerCase) != nrow(object@stand))
                        stop("length of 'lowerCase' != nrow of standardizing matrix")
                L2Fam <- eval(object@CallL2Fam)
                if(!identical(dim(L2Fam@param@trafo), dim(object@stand)))
                    stop("dimension of 'trafo' of 'param' != dimension of 'stand'")
                
                return(TRUE)
            })
# square integrable, conditionally centered (partial) IC 
# of contamination type
setClass("Av2CondContIC", 
            representation(clip = "numeric",
                           cent = "numeric",
                           stand = "numeric",
                           lowerCase = "OptionalNumeric",
                           neighborRadius = "numeric"), 
            prototype(name = "conditionally centered IC for average square conditional contamination neighborhoods",
                      Curve = EuclRandVarList(RealRandVariable(Map = list(function(x){x[1]*x[2]}), 
                                                    Domain = EuclideanSpace(dimension = 2))),
                      Risks = list(),
                      Infos = matrix(c(character(0),character(0)), ncol=2,
                                dimnames=list(character(0), c("method", "message"))),
                      CallL2Fam = call("L2RegTypeFamily"),
                      clip = Inf, 
                      cent = 0,
                      stand = 1,
                      lowerCase = NULL,
                      neighborRadius = 0),
            contains = "CondIC",
            validity = function(object){
                if(any(object@neighborRadius < 0)) # radius vector?!
                    stop("'neighborRadius' has to be in [0, Inf]")
                if(length(object@cent) != 1)
                    stop("length of 'cent' has to be 1")
                if(length(object@stand) != 1)
                    stop("length of 'stand' has to be 1")
                if(length(object@clip) != 1)
                    stop("length of clipping bound has to be 1")
                if(!is.null(object@lowerCase))
                    if(length(object@lowerCase) != 1)
                        stop("length of 'lowerCase' has to be 1")
                
                return(TRUE)
            })
# square integrable, conditionally centered (partial) IC 
# of total variation type
setClass("CondTotalVarIC", 
            representation(clipUp = "RealRandVariable",
                           clipLo = "RealRandVariable",
                           stand = "matrix",
                           lowerCase = "OptionalNumeric",
                           neighborRadius = "numeric",
                           neighborRadiusCurve = "function"), 
            prototype(name = "conditionally centered IC for average conditional contamination neighborhoods",
                      Curve = EuclRandVarList(RealRandVariable(Map = list(function(x){x[1]*x[2]}), 
                                               Domain = EuclideanSpace(dimension = 2))),
                      Risks = list(),
                      Infos = matrix(c(character(0),character(0)), ncol=2,
                                dimnames=list(character(0), c("method", "message"))),
                      CallL2Fam = call("L2RegTypeFamily"),
                      clipUp = RealRandVariable(Map = list(function(x){ Inf }),
                                                    Domain = EuclideanSpace(dimension = 1)),
                      clipLo = RealRandVariable(Map = list(function(x){ -Inf }),
                                                    Domain = EuclideanSpace(dimension = 1)),
                      stand = as.matrix(1),
                      lowerCase = NULL,
                      neighborRadius = 0,
                      neighborRadiusCurve = function(x){1}),
            contains = "CondIC",
            validity = function(object){
                if(any(object@neighborRadius < 0)) # radius vector?!
                    stop("'neighborRadius' has to be in [0, Inf]")
                if(length(formals(object@neighborRadiusCurve)) != 1)
                    stop("'neighborRadiusCurve' has to be a function of one argument")
                if(names(formals(object@neighborRadiusCurve)) != "x")
                    stop("'neighborRadiusCurve' has to be a function with argument name = 'x'")
                if(dimension(object@clipLo) != 1)
                    stop("dimension of lower clipping function has to be 1")
                if(dimension(object@clipUp) != 1)
                    stop("dimension of upper clipping function has to be 1")
                L2Fam <- eval(object@CallL2Fam)
                if(!identical(dim(L2Fam@param@trafo), dim(object@stand)))
                    stop("dimension of 'trafo' of 'param' != dimension of 'stand'")
                
                return(TRUE)
            })
# square integrable, conditionally centered (partial) IC 
# of total variation type
setClass("Av1CondTotalVarIC", 
            representation(clipUp = "numeric",
                           clipLo = "RealRandVariable",
                           stand = "matrix",
                           lowerCase = "OptionalNumeric",
                           neighborRadius = "numeric"), 
            prototype(name = "conditionally centered IC for average conditional contamination neighborhoods",
                      Curve = EuclRandVarList(RealRandVariable(Map = list(function(x){x[1]*x[2]}), 
                                               Domain = EuclideanSpace(dimension = 2))),
                      Risks = list(),
                      Infos = matrix(c(character(0),character(0)), ncol=2,
                                dimnames=list(character(0), c("method", "message"))),
                      CallL2Fam = call("L2RegTypeFamily"),
                      clipUp = Inf, 
                      clipLo = RealRandVariable(Map = list(function(x){ -Inf }),
                                                    Domain = EuclideanSpace(dimension = 1)),
                      stand = as.matrix(1),
                      lowerCase = NULL,
                      neighborRadius = 0),
            contains = "CondIC",
            validity = function(object){
                if(any(object@neighborRadius < 0)) # radius vector?!
                    stop("'neighborRadius' has to be in [0, Inf]")
                if(dimension(object@clipLo) != 1)
                    stop("dimension of lower clipping function has to be 1")
                if(length(object@clipUp) != 1)
                    stop("length of upper clipping bound has to be 1")
                L2Fam <- eval(object@CallL2Fam)
                if(!identical(dim(L2Fam@param@trafo), dim(object@stand)))
                    stop("dimension of 'trafo' of 'param' != dimension of 'stand'")
                
                return(TRUE)
            })
