## Generating function
Av2CondContIC <- function(name, CallL2Fam = call("L2RegTypeFamily"),
                   Curve = EuclRandVarList(RealRandVariable(Map = list(function(x){x[1]*x[2]}), 
                                                    Domain = EuclideanSpace(dimension = 2))), 
                   Risks, Infos, clip = Inf, stand = 1, cent = 0, 
                   lowerCase = NULL, neighborRadius = 0){
    if(missing(name))
        name <- "conditionally centered IC for average square conditional contamination neighborhoods"
    if(missing(Risks))
        Risks <- list()
    if(missing(Infos))
        Infos <- matrix(c(character(0),character(0)), ncol=2,
                    dimnames=list(character(0), c("method", "message")))
    
    return(new("Av2CondContIC", name = name, Curve = Curve, Risks = Risks, Infos = Infos,
               CallL2Fam = CallL2Fam, clip = clip, cent = cent, stand = stand, 
               lowerCase = lowerCase, neighborRadius = neighborRadius))
}

## generate IC
## for internal use only!
setMethod("generateIC", signature(neighbor = "Av2CondContNeighborhood", 
                                  L2Fam = "L2RegTypeFamily"),
    function(neighbor, L2Fam, res){
        A <- res$A
        z <- res$z
        b <- res$b
        d <- res$d
        ICfct <- vector(mode = "list", length = 1)
        L2 <- L2Fam@ErrorL2deriv[[1]]
        k <- dimension(img(L2Fam@RegDistr))
        K.inv <- solve(E(L2Fam@RegDistr, fun = function(x){ x %*% t(x) }))
        trafo <- L2Fam@param@trafo

        if(!is.null(d)){
            b0 <- b/sqrt(sum(diag(trafo %*% K.inv %*% t(trafo))))
            ICfct[[1]] <- function(x){ 
                                ind <- (L2(x[k+1]) != z) 
                                D <- matrix(D.vec, ncol = k)
                                K.inv <- matrix(K.vec, ncol = k)
                                b0*D %*% K.inv %*% x[1:k]*(sign(L2(x[k+1]) - z) + (1-ind)*d)
                          }
            body(ICfct[[1]]) <- substitute(
                                    { ind <- (L2(x[k+1]) != z) 
                                      D <- matrix(D.vec, ncol = k)
                                      K.inv <- matrix(K.vec, ncol = k)
                                      b0*D %*% K.inv %*% x[1:k]*(sign(L2(x[k+1]) - z) + (1-ind)*d) },
                                    list(L2 = L2@Map[[1]], D.vec = as.vector(trafo), 
                                         K.vec = as.vector(K.inv), b0 = b0, d = d, k = k))
        }else{
            c0 <- b/(A*sqrt(sum(diag(K.inv))))
            ICfct[[1]] <- function(x){ D <- matrix(D.vec, ncol = k)
                                       K.inv <- matrix(K.vec, ncol = k)
                                       A*D %*% K.inv %*% x[1:k]*(L2(x[k+1]) - z)*pmin(1, c0/abs(L2(x[k+1]) - z)) 
                          }
            body(ICfct[[1]]) <- substitute({ D <- matrix(D.vec, ncol = k)
                                             K.inv <- matrix(K.vec, ncol = k)
                                             A*D %*% K.inv %*% x[1:k]*(L2(x[k+1]) - z)*pmin(1, c0/abs(L2(x[k+1]) - z)) },
                                           list(L2 = L2@Map[[1]], D.vec = as.vector(trafo), 
                                                K.vec = as.vector(K.inv), c0 = c0, k = k))
        }

        return(Av2CondContIC(
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
                                         dimension = trunc(nrow(trafo)))),
                clip = b,
                cent = z,
                stand = A,
                lowerCase = d,
                neighborRadius = neighbor@radius,
                Risks = res$risk,
                Infos = matrix(res$info, ncol = 2, 
                            dimnames = list(character(0), c("method", "message")))))
    })

## Access methods
setMethod("clip", "Av2CondContIC", function(object) object@clip)
setMethod("cent", "Av2CondContIC", function(object) object@cent)
setMethod("stand", "Av2CondContIC", function(object) object@stand)
setMethod("lowerCase", "Av2CondContIC", function(object) object@lowerCase)
setMethod("neighborRadius", "Av2CondContIC", function(object) object@neighborRadius)

## replace methods
setReplaceMethod("clip", "Av2CondContIC", 
    function(object, value){ 
        stopifnot(is.numeric(value))
        L2Fam <- eval(object@CallL2Fam)
        res <- list(A = object@stand, z = object@cent, b = value, d = object@lowerCase,
                    risk = object@Risks, info = object@Infos)
        object <- generateIC(neighbor = Av2CondContNeighborhood(radius = object@neighborRadius), 
                             L2Fam = L2Fam, res = res)
        addInfo(object) <- c("clip<-", "The clipping bound has been changed")
        addInfo(object) <- c("clip<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
setReplaceMethod("cent", "Av2CondContIC", 
    function(object, value){ 
        stopifnot(is.numeric(value))
        L2Fam <- eval(object@CallL2Fam)
        res <- list(A = object@stand, z = value, b = object@clip, d = object@lowerCase,
                    risk = object@Risks, info = object@Infos)
        object <- generateIC(neighbor = Av2CondContNeighborhood(radius = object@neighborRadius), 
                             L2Fam = L2Fam, res = res)
        addInfo(object) <- c("cent<-", "The centering constant has been changed")
        addInfo(object) <- c("cent<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
setReplaceMethod("stand", "Av2CondContIC", 
    function(object, value){ 
        stopifnot(is.numeric(value))
        L2Fam <- eval(object@CallL2Fam)
        res <- list(A = value, z = object@cent, b = object@clip, d = object@lowerCase,
                    risk = object@Risks, info = object@Infos)
        object <- generateIC(neighbor = Av2CondContNeighborhood(radius = object@neighborRadius), 
                             L2Fam = L2Fam, res = res)
        addInfo(object) <- c("stand<-", "The standardizing matrix has been changed")
        addInfo(object) <- c("stand<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
setReplaceMethod("lowerCase", "Av2CondContIC", 
    function(object, value){ 
        stopifnot(is.null(value)||is.numeric(value))
        L2Fam <- eval(object@CallL2Fam)
        res <- list(A = object@stand, z = object@cent, b = object@clip, d = value,
                    risk = object@Risks, info = object@Infos)
        object <- generateIC(neighbor = Av2CondContNeighborhood(radius = object@neighborRadius), 
                             L2Fam = L2Fam, res = res)
        addInfo(object) <- c("lowerCase<-", "The slot 'lowerCase' has been changed")
        addInfo(object) <- c("lowerCase<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
setReplaceMethod("neighborRadius", "Av2CondContIC", 
    function(object, value){ 
        object@neighborRadius <- value
        if(any(value < 0)) # radius vector?!
            stop("'value' has to be in [0, Inf]")
        addInfo(object) <- c("neighborRadius<-", "The slot 'neighborRadius' has been changed")
        addInfo(object) <- c("neighborRadius<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
setReplaceMethod("CallL2Fam", "Av2CondContIC",
    function(object, value){ 
        L2Fam <- eval(value)
        res <- list(A = object@stand, z = object@cent, b = object@clip, d = object@lowerCase,
                    risk = object@Risks, info = object@Infos)
        object <- generateIC(neighbor = Av2CondContNeighborhood(radius = object@neighborRadius), 
                             L2Fam = L2Fam, res = res)
        addInfo(object) <- c("CallL2Fam<-", "The slot 'CallL2Fam' has been changed")
        addInfo(object) <- c("CallL2Fam<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
