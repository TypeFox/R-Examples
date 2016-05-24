## Generating function
Av1CondTotalVarIC <- function(name, CallL2Fam = call("L2RegTypeFamily"),
                   Curve = EuclRandVarList(RealRandVariable(Map = list(function(x){x[1]*x[2]}), 
                                                    Domain = EuclideanSpace(dimension = 2))), 
                   Risks, Infos, clipUp = Inf, stand = as.matrix(1), 
                   clipLo = RealRandVariable(Map = list(function(x){ -Inf }),
                                             Domain = EuclideanSpace(dimension = 1)), 
                   lowerCase = NULL, neighborRadius = 0){
    if(missing(name))
        name <- "conditionally centered IC for average conditional total variation neighborhoods"
    if(missing(Risks))
        Risks <- list()
    if(missing(Infos))
        Infos <- matrix(c(character(0),character(0)), ncol=2,
                    dimnames=list(character(0), c("method", "message")))
    
    return(new("Av1CondTotalVarIC", name = name, Curve = Curve, Risks = Risks, Infos = Infos,
               CallL2Fam = CallL2Fam, clipUp = clipUp, clipLo = clipLo, stand = stand, 
               lowerCase = lowerCase, neighborRadius = neighborRadius))
}

## generate IC
## for internal use only!
setMethod("generateIC", signature(neighbor = "Av1CondTotalVarNeighborhood", 
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
            ICfct[[1]] <- function(x){ ind1 <- (L2(x[k+1]) > 0); ind2 <- (L2(x[k+1]) < 0)
                                       A <- matrix(A.vec, ncol = k)
                                       Y <- as.vector(A %*% x[1:k])
                                       v <- sqrt(sum(Y^2))
                                       ax <- a(x[1:k])
                                       Y/v*((ax+b)*ind1 + ax*ind2)
                          }
            body(ICfct[[1]]) <- substitute({ ind1 <- (L2(x[k+1]) > 0); ind2 <- (L2(x[k+1]) < 0)
                                             A <- matrix(A.vec, ncol = k)
                                             Y <- as.vector(A %*% x[1:k])
                                             v <- sqrt(sum(Y^2))
                                             ax <- a(x[1:k])
                                             Y/v*((ax+b)*ind1 + ax*ind2) },
                                           list(A.vec = as.vector(A), L2 = L2@Map[[1]], a = a@Map[[1]], 
                                                b = b, k = k))
        }else{
            if(b == Inf){
                ICfct[[1]]<- function(x){ A <- matrix(A.vec, ncol = k)
                                          v <- as.vector(sqrt(sum((A %*% x[1:k])^2)))
                                          ax <- a(x[1:k])
                                          if(ax == -Inf) 
                                              as.vector(A %*% x[1:k])*L2(x[k+1])
                                          else 
                                              as.vector(A %*% x[1:k])*max(a(x[1:k])/v, L2(x[k+1])) }
                body(ICfct[[1]]) <- substitute({ A <- matrix(A.vec, ncol = k)
                                                 v <- as.vector(sqrt(sum((A %*% x[1:k])^2)))
                                                 ax <- a(x[1:k])
                                                 if(ax == -Inf) 
                                                     as.vector(A %*% x[1:k])*L2(x[k+1])
                                                 else 
                                                     as.vector(A %*% x[1:k])*max(a(x[1:k])/v, L2(x[k+1])) }, 
                                               list(A.vec = as.vector(A), L2 = L2@Map[[1]], a = a@Map[[1]], 
                                                    b = b, k = k))
            }else{
                ICfct[[1]] <- function(x){ A <- matrix(A.vec, ncol = k)
                                           v <- as.vector(sqrt(sum((A %*% x[1:k])^2)))
                                           ax <- a(x[1:k])
                                           as.vector(A %*% x[1:k])*min(max(a(x[1:k])/v, L2(x[k+1])), (ax+b)/v) }
                body(ICfct[[1]]) <- substitute({ A <- matrix(A.vec, ncol = k)
                                                 v <- as.vector(sqrt(sum((A %*% x[1:k])^2)))
                                                 ax <- a(x[1:k])
                                                 as.vector(A %*% x[1:k])*min(max(ax/v, L2(x[k+1])), (ax+b)/v) },
                                               list(A.vec = as.vector(A), L2 = L2@Map[[1]], a = a@Map[[1]], 
                                                    b = b, k = k))
            }
        }
        return(Av1CondTotalVarIC(
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
                Risks = res$risk,
                Infos = matrix(res$info, ncol = 2, 
                            dimnames = list(character(0), c("method", "message")))))
    })

## Access methods
setMethod("clipUp", "Av1CondTotalVarIC", function(object) object@clipUp)
setMethod("clipLo", "Av1CondTotalVarIC", function(object) object@clipLo)
setMethod("stand", "Av1CondTotalVarIC", function(object) object@stand)
setMethod("lowerCase", "Av1CondTotalVarIC", function(object) object@lowerCase)
setMethod("neighborRadius", "Av1CondTotalVarIC", function(object) object@neighborRadius)

## replace methods
setReplaceMethod("clipUp", "Av1CondTotalVarIC", 
    function(object, value){ 
        stopifnot(is.numeric(value))
        L2Fam <- eval(object@CallL2Fam)
        res <- list(A = object@stand, a = object@clipLo, b = value, d = object@lowerCase,
                    risk = object@Risks, info = object@Infos)
        object <- generateIC(neighbor = Av1CondTotalVarNeighborhood(radius = object@neighborRadius), 
                             L2Fam = L2Fam, res = res)
        addInfo(object) <- c("clipUp<-", "The clipping bound has been changed")
        addInfo(object) <- c("clipUp<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
setReplaceMethod("clipLo", "Av1CondTotalVarIC", 
    function(object, value){ 
        stopifnot(is(value, "RealRandVariable"))
        L2Fam <- eval(object@CallL2Fam)
        res <- list(A = object@stand, a = value, b = object@clipUp, d = object@lowerCase,
                    risk = object@Risks, info = object@Infos)
        object <- generateIC(neighbor = Av1CondTotalVarNeighborhood(radius = object@neighborRadius), 
                             L2Fam = L2Fam, res = res)
        addInfo(object) <- c("clipLo<-", "The centering constant has been changed")
        addInfo(object) <- c("clipLo<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
setReplaceMethod("stand", "Av1CondTotalVarIC", 
    function(object, value){ 
        stopifnot(is.matrix(value))
        L2Fam <- eval(object@CallL2Fam)
        res <- list(A = value, a = object@clipLo, b = object@clipUp, d = object@lowerCase,
                    risk = object@Risks, info = object@Infos)
        object <- generateIC(neighbor = Av1CondTotalVarNeighborhood(radius = object@neighborRadius), 
                             L2Fam = L2Fam, res = res)
        addInfo(object) <- c("stand<-", "The standardizing matrix has been changed")
        addInfo(object) <- c("stand<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
setReplaceMethod("lowerCase", "Av1CondTotalVarIC", 
    function(object, value){ 
        stopifnot(is.null(value)||is.numeric(value))
        L2Fam <- eval(object@CallL2Fam)
        res <- list(A = object@stand, a = object@clipLo, b = object@clipUp, d = value,
                    risk = object@Risks, info = object@Infos)
        object <- generateIC(neighbor = Av1CondTotalVarNeighborhood(radius = object@neighborRadius), 
                             L2Fam = L2Fam, res = res)
        addInfo(object) <- c("lowerCase<-", "The slot 'lowerCase' has been changed")
        addInfo(object) <- c("lowerCase<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
setReplaceMethod("neighborRadius", "Av1CondTotalVarIC", 
    function(object, value){ 
        object@neighborRadius <- value
        if(any(value < 0)) # radius vector?!
            stop("'value' has to be in [0, Inf]")
        addInfo(object) <- c("neighborRadius<-", "The slot 'neighborRadius' has been changed")
        addInfo(object) <- c("neighborRadius<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
setReplaceMethod("CallL2Fam", "Av1CondTotalVarIC",
    function(object, value){ 
        L2Fam <- eval(value)
        res <- list(A = object@stand, a = object@clipLo, b = object@clipUp, d = object@lowerCase,
                    risk = object@Risks, info = object@Infos)
        object <- generateIC(neighbor = Av1CondTotalVarNeighborhood(radius = object@neighborRadius), 
                             L2Fam = L2Fam, res = res)
        addInfo(object) <- c("CallL2Fam<-", "The slot 'CallL2Fam' has been changed")
        addInfo(object) <- c("CallL2Fam<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
