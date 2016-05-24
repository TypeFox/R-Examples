## Generating function
CondContIC <- function(name, CallL2Fam = call("L2RegTypeFamily"),
                   Curve = EuclRandVarList(RealRandVariable(Map = list(function(x){x[1]*x[2]}), 
                                                    Domain = EuclideanSpace(dimension = 2))), 
                   Risks, Infos, 
                   clip = RealRandVariable(Map = list(function(x){ Inf }), Domain = Reals()), 
                   stand = as.matrix(1), 
                   cent = EuclRandVarList(RealRandVariable(Map = list(function(x){numeric(length(x))}),
                                                    Domain = EuclideanSpace(dimension = 2))), 
                   lowerCase = NULL, neighborRadius = 0, neighborRadiusCurve = function(x){1}){
    if(missing(name))
        name <- "conditionally centered IC for average conditional contamination neighborhoods"
    if(missing(Risks))
        Risks <- list()
    if(missing(Infos))
        Infos <- matrix(c(character(0),character(0)), ncol=2,
                    dimnames=list(character(0), c("method", "message")))
    
    return(new("CondContIC", name = name, Curve = Curve, Risks = Risks, Infos = Infos,
               CallL2Fam = CallL2Fam, clip = clip, cent = cent, stand = stand, 
               lowerCase = lowerCase, neighborRadius = neighborRadius, 
               neighborRadiusCurve = neighborRadiusCurve))
}

## generate IC
## for internal use only!
setMethod("generateIC", signature(neighbor = "CondContNeighborhood", 
                                  L2Fam = "L2RegTypeFamily"),
    function(neighbor, L2Fam, res){
        A <- res$A
        a <- res$a
        b <- res$b
        d <- res$d
        nrvalues <- nrow(A)
        ICfct <- vector(mode = "list", length = nrvalues)
        a1 <- as(diag(nrvalues) %*% a, "EuclRandVariable")
        Y <- as(A %*% L2Fam@L2deriv, "EuclRandVariable") - a1
        k <- dimension(img(L2Fam@RegDistr))
        if(nrvalues == 1){
            if(!is.null(d)){
                ICfct[[1]] <- function(x){ 
                                    ind <- (Y(x) != 0) 
                                    b(x[1:k])*(ind*Y(x)/(ind*absY(x) + (1-ind)) + zi*(1-ind)*d)
                              }
                body(ICfct[[1]]) <- substitute(
                                        { ind <- (Y(x) != 0) 
                                          b(x[1:k])*(ind*Y(x)/(ind*absY(x) + (1-ind)) + zi*(1-ind)*d) },
                                        list(Y = Y@Map[[1]], absY = abs(Y)@Map[[1]], b = b@Map[[1]], d = d, 
                                             zi = sign(L2Fam@param@trafo), k = k))
            }else{
                ICfct[[1]] <- function(x){ Y(x)*pmin(1, b(x[1:k])/absY(x)) }
                body(ICfct[[1]]) <- substitute({ Y(x)*pmin(1, b(x[1:k])/absY(x)) },
                                                 list(Y = Y@Map[[1]], absY = abs(Y)@Map[[1]], 
                                                      b = b@Map[[1]], k = k))
            }
        }else{
            absY <- sqrt(Y %*% Y)
            if(!is.null(d))
                for(i in 1:nrvalues){
                    ICfct[[i]] <- function(x){ ind <- (Yi(x) != 0) ; ind*b(x[1:k])*Yi(x)/absY(x) + (1-ind)*d }
                    body(ICfct[[i]]) <- substitute({ ind <- (Yi(x) != 0) ; ind*b(x[1:k])*Yi(x)/absY(x) + (1-ind)*d },
                                                 list(Yi = Y@Map[[i]], absY = absY@Map[[1]], b = b@Map[[1]], 
                                                      d = d[i], k = k))
                }
            else
                for(i in 1:nrvalues){
                    ICfct[[i]] <- function(x){ Yi(x)*pmin(1, b(x[1:k])/absY(x)) }
                    body(ICfct[[i]]) <- substitute({ Yi(x)*pmin(1, b(x[1:k])/absY(x)) },
                                                 list(Yi = Y@Map[[i]], absY = absY@Map[[1]], b = b@Map[[1]], k = k))
                }
        }
        return(CondContIC(
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
                Curve = EuclRandVarList(EuclRandVariable(Map = ICfct, Domain = Y@Domain, 
                                         Range = Y@Range)),
                clip = b,
                cent = a,
                stand = A,
                lowerCase = d,
                neighborRadius = neighbor@radius,
                neighborRadiusCurve = neighbor@radiusCurve,
                Risks = res$risk,
                Infos = matrix(res$info, ncol = 2, 
                            dimnames = list(character(0), c("method", "message")))))
    })

## Access methods
setMethod("clip", "CondContIC", function(object) object@clip)
setMethod("cent", "CondContIC", function(object) object@cent)
setMethod("stand", "CondContIC", function(object) object@stand)
setMethod("lowerCase", "CondContIC", function(object) object@lowerCase)
setMethod("neighborRadius", "CondContIC", function(object) object@neighborRadius)
setMethod("neighborRadiusCurve", "CondContIC", function(object) object@neighborRadiusCurve)

## replace methods
setReplaceMethod("clip", "CondContIC", 
    function(object, value){ 
        stopifnot(is(value, "RealRandVariable"))
        L2Fam <- eval(object@CallL2Fam)
        res <- list(A = object@stand, a = object@cent, b = value, d = object@lowerCase,
                    risk = object@Risks, info = object@Infos)
        object <- generateIC(neighbor = CondContNeighborhood(radius = object@neighborRadius, 
                                                        radiusCurve = object@neighborRadiusCurve), 
                             L2Fam = L2Fam, res = res)
        addInfo(object) <- c("clip<-", "The clipping bound has been changed")
        addInfo(object) <- c("clip<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
setReplaceMethod("cent", "CondContIC", 
    function(object, value){ 
        stopifnot(is(value, "EuclRandVarList"))
        L2Fam <- eval(object@CallL2Fam)
        res <- list(A = object@stand, a = value, b = object@clip, d = object@lowerCase,
                    risk = object@Risks, info = object@Infos)
        object <- generateIC(neighbor = CondContNeighborhood(radius = object@neighborRadius, 
                                                        radiusCurve = object@neighborRadiusCurve), 
                             L2Fam = L2Fam, res = res)
        addInfo(object) <- c("cent<-", "The centering constant has been changed")
        addInfo(object) <- c("cent<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
setReplaceMethod("stand", "CondContIC", 
    function(object, value){ 
        stopifnot(is.matrix(value))
        L2Fam <- eval(object@CallL2Fam)
        res <- list(A = value, a = object@cent, b = object@clip, d = object@lowerCase,
                    risk = object@Risks, info = object@Infos)
        object <- generateIC(neighbor = CondContNeighborhood(radius = object@neighborRadius, 
                                                        radiusCurve = object@neighborRadiusCurve), 
                             L2Fam = L2Fam, res = res)
        addInfo(object) <- c("stand<-", "The standardizing matrix has been changed")
        addInfo(object) <- c("stand<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
setReplaceMethod("lowerCase", "CondContIC", 
    function(object, value){ 
        stopifnot(is.null(value)||is.numeric(value))
        L2Fam <- eval(object@CallL2Fam)
        res <- list(A = object@stand, a = object@cent, b = object@clip, d = value,
                    risk = object@Risks, info = object@Infos)
        object <- generateIC(neighbor = CondContNeighborhood(radius = object@neighborRadius, 
                                                        radiusCurve = object@neighborRadiusCurve), 
                             L2Fam = L2Fam, res = res)
        addInfo(object) <- c("lowerCase<-", "The slot 'lowerCase' has been changed")
        addInfo(object) <- c("lowerCase<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
setReplaceMethod("neighborRadius", "CondContIC", 
    function(object, value){ 
        object@neighborRadius <- value
        if(any(value < 0)) # radius vector?!
            stop("'value' has to be in [0, Inf]")
        addInfo(object) <- c("neighborRadius<-", "The slot 'neighborRadius' has been changed")
        addInfo(object) <- c("neighborRadius<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
setReplaceMethod("neighborRadiusCurve", "CondContIC", 
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
setReplaceMethod("CallL2Fam", "CondContIC",
    function(object, value){ 
        L2Fam <- eval(value)
        res <- list(A = object@stand, a = object@cent, b = object@clip, d = object@lowerCase,
                    risk = object@Risks, info = object@Infos)
        object <- generateIC(neighbor = CondContNeighborhood(radius = object@neighborRadius, 
                                                        radiusCurve = object@neighborRadiusCurve), 
                             L2Fam = L2Fam, res = res)
        addInfo(object) <- c("CallL2Fam<-", "The slot 'CallL2Fam' has been changed")
        addInfo(object) <- c("CallL2Fam<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
