## Generating function
CondTotalVarIC <- function(name, CallL2Fam = call("L2RegTypeFamily"),
                   Curve = EuclRandVarList(RealRandVariable(Map = list(function(x){x[1]*x[2]}), 
                                                    Domain = EuclideanSpace(dimension = 2))), 
                   Risks, Infos, 
                   clipUp = RealRandVariable(Map = list(function(x){ Inf }), Domain = Reals()),
                   stand = as.matrix(1), 
                   clipLo = RealRandVariable(Map = list(function(x){ -Inf }), Domain = Reals()), 
                   lowerCase = NULL, neighborRadius = 0, neighborRadiusCurve = function(x){1}){
    if(missing(name))
        name <- "conditionally centered IC for average conditional total variation neighborhoods"
    if(missing(Risks))
        Risks <- list()
    if(missing(Infos))
        Infos <- matrix(c(character(0),character(0)), ncol=2,
                    dimnames=list(character(0), c("method", "message")))
    
    return(new("CondTotalVarIC", name = name, Curve = Curve, Risks = Risks, Infos = Infos,
               CallL2Fam = CallL2Fam, clipUp = clipUp, clipLo = clipLo, stand = stand, 
               lowerCase = lowerCase, neighborRadius = neighborRadius, 
               neighborRadiusCurve = neighborRadiusCurve))
}

## generate IC
## for internal use only!
setMethod("generateIC", signature(neighbor = "CondTotalVarNeighborhood", 
                                  L2Fam = "L2RegTypeFamily"),
    function(neighbor, L2Fam, res){
        A <- res$A
        a <- res$a
        b <- res$b
        d <- res$d
        ICfct <- vector(mode = "list", length = 1)
        L2 <- L2Fam@ErrorL2deriv[[1]]
        k <- dimension(img(L2Fam@RegDistr))
        if(!is.null(d)){
            ICfct[[1]] <- function(x){ A <- matrix(A.vec, ncol = k)
                                       Y <- as.vector(A %*% x[1:k]) * L2(x[k+1])
                                       ind1 <- (Y > 0); ind2 <- (Y < 0)
                                       b(x[1:k])*ind1 + a(x[1:k])*ind2
                          }
            body(ICfct[[1]]) <- substitute({ A <- matrix(A.vec, ncol = k)
                                             Y <- as.vector(A %*% x[1:k]) * L2(x[k+1])
                                             ind1 <- (Y > 0); ind2 <- (Y < 0)
                                             b(x[1:k])*ind1 + a(x[1:k])*ind2 },
                                           list(A.vec = as.vector(A), L2 = L2@Map[[1]], a = a@Map[[1]], 
                                                b = b@Map[[1]], k = k))
        }else{
            ICfct[[1]]<- function(x){ A <- matrix(A.vec, ncol = k)
                                      Y <- as.vector(A %*% x[1:k]) * L2(x[k+1])
                                      min(max(a(x[1:k]), Y), b(x[1:k])) }
            body(ICfct[[1]]) <- substitute({ A <- matrix(A.vec, ncol = k)
                                             Y <- as.vector(A %*% x[1:k]) * L2(x[k+1])
                                             min(max(a(x[1:k]), Y), b(x[1:k])) }, 
                                           list(A.vec = as.vector(A), L2 = L2@Map[[1]], a = a@Map[[1]], 
                                                b = b@Map[[1]], k = k))
        }
        return(CondTotalVarIC(
                name = "conditionally centered IC of contamination type", 
                CallL2Fam = call("L2RegTypeFamily", 
                                name = L2Fam@name,
                                distribution = L2Fam@distribution,  
                                param = L2Fam@param,
                                props = L2Fam@props,
                                L2deriv = L2Fam@L2deriv,
                                ErrorDistr = L2Fam@ErrorDistr,
                                ErrorSymm = L2Fam@ErrorSymm,
                                RegDistr = L2Fam@RegDistr,
                                RegSymm = L2Fam@RegSymm,
                                Regressor = L2Fam@Regressor,
                                ErrorL2deriv = L2Fam@ErrorL2deriv,
                                ErrorL2derivDistr = L2Fam@ErrorL2derivDistr,
                                ErrorL2derivSymm = L2Fam@ErrorL2derivSymm,
                                FisherInfo = L2Fam@FisherInfo),
                Curve = EuclRandVarList(EuclRandVariable(Map = ICfct, Domain = L2Fam@L2deriv[[1]]@Domain, 
                                        dimension = trunc(nrow(L2Fam@param@trafo)))),
                clipUp = b,
                clipLo = a,
                stand = A,
                lowerCase = d,
                neighborRadius = neighbor@radius,
                neighborRadiusCurve = neighbor@radiusCurve,
                Risks = res$risk,
                Infos = matrix(res$info, ncol = 2, 
                            dimnames = list(character(0), c("method", "message")))))
    })

## Access methods
setMethod("clipUp", "CondTotalVarIC", function(object) object@clipUp)
setMethod("clipLo", "CondTotalVarIC", function(object) object@clipLo)
setMethod("stand", "CondTotalVarIC", function(object) object@stand)
setMethod("lowerCase", "CondTotalVarIC", function(object) object@lowerCase)
setMethod("neighborRadius", "CondTotalVarIC", function(object) object@neighborRadius)
setMethod("neighborRadiusCurve", "CondTotalVarIC", function(object) object@neighborRadiusCurve)

## replace methods
setReplaceMethod("clipUp", "CondTotalVarIC", 
    function(object, value){ 
        stopifnot(is(value, "RealRandVariable"))
        L2Fam <- eval(object@CallL2Fam)
        res <- list(A = object@stand, a = object@clipLo, b = value, d = object@lowerCase,
                    risk = object@Risks, info = object@Infos)
        object <- generateIC(neighbor = CondTotalVarNeighborhood(radius = object@neighborRadius,
                                                            radiusCurve = object@neighborRadiusCurve), 
                             L2Fam = L2Fam, res = res)
        addInfo(object) <- c("clipUp<-", "The clipping bound has been changed")
        addInfo(object) <- c("clipUp<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
setReplaceMethod("clipLo", "CondTotalVarIC", 
    function(object, value){ 
        stopifnot(is(value, "RealRandVariable"))
        L2Fam <- eval(object@CallL2Fam)
        res <- list(A = object@stand, a = value, b = object@clipUp, d = object@lowerCase,
                    risk = object@Risks, info = object@Infos)
        object <- generateIC(neighbor = CondTotalVarNeighborhood(radius = object@neighborRadius,
                                                            radiusCurve = object@neighborRadiusCurve), 
                             L2Fam = L2Fam, res = res)
        addInfo(object) <- c("clipLo<-", "The centering constant has been changed")
        addInfo(object) <- c("clipLo<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
setReplaceMethod("stand", "CondTotalVarIC", 
    function(object, value){ 
        stopifnot(is.matrix(value))
        L2Fam <- eval(object@CallL2Fam)
        res <- list(A = value, a = object@clipLo, b = object@clipUp, d = object@lowerCase,
                    risk = object@Risks, info = object@Infos)
        object <- generateIC(neighbor = CondTotalVarNeighborhood(radius = object@neighborRadius, 
                                                            radiusCurve = object@neighborRadiusCurve), 
                             L2Fam = L2Fam, res = res)
        addInfo(object) <- c("stand<-", "The standardizing matrix has been changed")
        addInfo(object) <- c("stand<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
setReplaceMethod("lowerCase", "CondTotalVarIC", 
    function(object, value){ 
        stopifnot(is.null(value)||is.numeric(value))
        L2Fam <- eval(object@CallL2Fam)
        res <- list(A = object@stand, a = object@clipLo, b = object@clipUp, d = value,
                    risk = object@Risks, info = object@Infos)
        object <- generateIC(neighbor = CondTotalVarNeighborhood(radius = object@neighborRadius, 
                                                            radiusCurve = object@neighborRadiusCurve), 
                             L2Fam = L2Fam, res = res)
        addInfo(object) <- c("lowerCase<-", "The slot 'lowerCase' has been changed")
        addInfo(object) <- c("lowerCase<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
setReplaceMethod("neighborRadius", "CondTotalVarIC", 
    function(object, value){ 
        object@neighborRadius <- value
        if(any(value < 0)) # radius vector?!
            stop("'value' has to be in [0, Inf]")
        addInfo(object) <- c("neighborRadius<-", "The slot 'neighborRadius' has been changed")
        addInfo(object) <- c("neighborRadius<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
setReplaceMethod("neighborRadiusCurve", "CondTotalVarIC", 
    function(object, value){ 
        object@neighborRadiusCurve <- value
        if(length(formals(value)) != 1)
            stop("'value' has to be a function of one argument")
        if(names(formals(value)) != "x")
            stop("'value' has to be a function with argument name = 'x'")
        addInfo(object) <- c("neighborRadiusCurve<-", "The slot 'neighborRadiusCurve' has been changed")
        addInfo(object) <- c("neighborRadiusCurve<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
setReplaceMethod("CallL2Fam", "CondTotalVarIC",
    function(object, value){ 
        L2Fam <- eval(value)
        res <- list(A = object@stand, a = object@clipLo, b = object@clipUp, d = object@lowerCase,
                    risk = object@Risks, info = object@Infos)
        object <- generateIC(neighbor = CondTotalVarNeighborhood(radius = object@neighborRadius, 
                                                            radiusCurve = object@neighborRadiusCurve), 
                             L2Fam = L2Fam, res = res)
        addInfo(object) <- c("CallL2Fam<-", "The slot 'CallL2Fam' has been changed")
        addInfo(object) <- c("CallL2Fam<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
