## Generating function
TotalVarIC <- function(name, CallL2Fam = call("L2ParamFamily"),
                    Curve = EuclRandVarList(RealRandVariable(Map = c(function(x){x}), Domain = Reals())), 
                    Risks, Infos, clipLo = -Inf, clipUp = Inf, stand = as.matrix(1),
                    lowerCase = NULL, neighborRadius = 0){
    if(missing(name))
        name <- "IC of total variation type"
    if(missing(Risks))
        Risks <- list()
    if(missing(Infos))
        Infos <- matrix(c(character(0),character(0)), ncol=2,
                    dimnames=list(character(0), c("method", "message")))

    if(any(neighborRadius < 0)) # radius vector?!
        stop("'neighborRadius' has to be in [0, Inf]")
    if((length(clipLo) != 1) && (length(clipLo) != length(Curve)))
        stop("length of lower clipping bound != 1 and != length of 'Curve'")
    L2Fam <- eval(CallL2Fam)
    if(!identical(dim(L2Fam@param@trafo), dim(stand)))
        stop(paste("dimension of 'trafo' of 'param' != dimension of 'stand'"))
    
    IC1 <- new("TotalVarIC")
    IC1@name <- name
    IC1@Curve <- Curve
    IC1@Risks <- Risks
    IC1@Infos <- Infos
    IC1@CallL2Fam <- CallL2Fam
    IC1@clipLo <- clipLo
    IC1@clipUp <- clipUp
    IC1@stand <- stand
    IC1@lowerCase <- lowerCase
    IC1@neighborRadius <- neighborRadius
    
    return(IC1)
}

## generate IC
## for internal use only!
setMethod("generateIC", signature(neighbor = "TotalVarNeighborhood", 
                                  L2Fam = "L2ParamFamily"),
    function(neighbor, L2Fam, res){
        A <- res$A
        a <- sign(as.vector(A))*res$a
        b <- res$b
        ICfct <- vector(mode = "list", length = 1)
        Y <- as(A %*% L2Fam@L2deriv, "EuclRandVariable")
        if(!is.null(res$d)){
            a <- as.vector(A)*a
            ICfct[[1]] <- function(x){ ind1 <- (Y(x) > 0); ind2 <- (Y(x) < 0)
                                       (a+b)*ind1 + a*ind2 }
            body(ICfct[[1]]) <- substitute({ ind1 <- (Y(x) > 0); ind2 <- (Y(x) < 0)
                                             (a+b)*ind1 + a*ind2 },
                                             list(Y = Y@Map[[1]], a = a, b = b))
        }else{
            if((a == -Inf)&&(b == Inf)){
                ICfct[[1]]<- function(x){ Y(x) }
                body(ICfct[[1]]) <- substitute({ Y(x) }, list(Y = Y@Map[[1]]))
            }else{
                ICfct[[1]] <- function(x){ min(max(a, Y(x)), a+b) }
                body(ICfct[[1]]) <- substitute({ min(max(a, Y(x)), a+b) },
                                                 list(Y = Y@Map[[1]], a = a, b = b))
            }
        }
        if((a == -Inf) & (b == Inf))
            clipUp <- Inf
        else
            clipUp <- a + b
        return(TotalVarIC(
                name = "IC of total variation type", 
                CallL2Fam = call("L2ParamFamily", 
                                name = L2Fam@name,
                                distribution = L2Fam@distribution, 
                                distrSymm = L2Fam@distrSymm, 
                                param = L2Fam@param,
                                props = L2Fam@props,
                                L2deriv = L2Fam@L2deriv,
                                L2derivSymm = L2Fam@L2derivSymm,
                                L2derivDistr = L2Fam@L2derivDistr,
                                L2derivDistrSymm = L2Fam@L2derivDistrSymm,
                                FisherInfo = L2Fam@FisherInfo),
                Curve = EuclRandVarList(EuclRandVariable(Map = ICfct, Domain = Y@Domain, 
                            Range = Y@Range)),
                clipUp = clipUp,
                clipLo = a,
                stand = A,
                lowerCase = res$d, 
                neighborRadius = neighbor@radius,
                Risks = res$risk,
                Infos = matrix(res$info, ncol = 2, 
                            dimnames = list(character(0), c("method", "message")))))
    })

## Access methods
setMethod("clipLo", "TotalVarIC", function(object) object@clipLo)
setMethod("clipUp", "TotalVarIC", function(object) object@clipUp)
setMethod("stand", "TotalVarIC", function(object) object@stand)
setMethod("lowerCase", "TotalVarIC", function(object) object@lowerCase)
setMethod("neighborRadius", "TotalVarIC", function(object) object@neighborRadius)

## Replace methods
setReplaceMethod("clipLo", "TotalVarIC", 
    function(object, value){ 
        stopifnot(is.numeric(value))
        L2Fam <- eval(object@CallL2Fam)
        res <- list(A = object@stand, a = value, b = object@clipUp-value, 
                    d = object@lowerCase, risk = object@Risks, info = object@Infos)
        object <- generateIC(neighbor = TotalVarNeighborhood(radius = object@neighborRadius), 
                             L2Fam = L2Fam, res = res)
        addInfo(object) <- c("clipLo<-", "The lower clipping bound has been changed")
        addInfo(object) <- c("clipLo<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
setReplaceMethod("clipUp", "TotalVarIC", 
    function(object, value){ 
        stopifnot(is.numeric(value))
        L2Fam <- eval(object@CallL2Fam)
        res <- list(A = object@stand, a = object@clipLo, b = value-object@clipLo, 
                    d = object@lowerCase, risk = object@Risks, info = object@Infos)
        object <- generateIC(neighbor = TotalVarNeighborhood(radius = object@neighborRadius), 
                             L2Fam = L2Fam, res = res)
        addInfo(object) <- c("clipUp<-", "The upper clipping bound has been changed")
        addInfo(object) <- c("clipUp<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
setReplaceMethod("stand", "TotalVarIC", 
    function(object, value){ 
        stopifnot(is.matrix(value))
        L2Fam <- eval(object@CallL2Fam)
        res <- list(A = value, a = object@clipLo, b = object@clipUp-object@clipLo, 
                    d = object@lowerCase, risk = object@Risks, info = object@Infos)
        object <- generateIC(neighbor = TotalVarNeighborhood(radius = object@neighborRadius), 
                             L2Fam = L2Fam, res = res)
        addInfo(object) <- c("stand<-", "The standardizing matrix has been changed")
        addInfo(object) <- c("stand<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
setReplaceMethod("lowerCase", "TotalVarIC", 
    function(object, value){ 
        stopifnot(is.null(value)||is.numeric(value))
        L2Fam <- eval(object@CallL2Fam)
        res <- list(A = object@stand, a = object@clipLo, b = object@clipUp-object@clipLo, 
                    d = value, risk = object@Risks, info = object@Infos)
        object <- generateIC(neighbor = TotalVarNeighborhood(radius = object@neighborRadius), 
                             L2Fam = L2Fam, res = res)
        addInfo(object) <- c("lowerCase<-", "The slot 'lowerCase' has been changed")
        addInfo(object) <- c("lowerCase<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
setReplaceMethod("neighborRadius", "TotalVarIC", 
    function(object, value){ 
        object@neighborRadius <- value
        if(any(value < 0)) # radius vector?!
            stop("'value' has to be in [0, Inf]")
        addInfo(object) <- c("neighborRadius<-", "The slot 'neighborRadius' has been changed")
        addInfo(object) <- c("neighborRadius<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
setReplaceMethod("CallL2Fam", "TotalVarIC",
    function(object, value){ 
        L2Fam <- eval(value)
        res <- list(A = object@stand, a = object@clipLo, b = object@clipUp-object@clipLo, 
                    d = object@lowerCase, risk = object@Risks, info = object@Infos)
        object <- generateIC(neighbor = TotalVarNeighborhood(radius = object@neighborRadius), 
                             L2Fam = L2Fam, res = res)
        addInfo(object) <- c("CallL2Fam<-", "The slot 'CallL2Fam' has been changed")
        addInfo(object) <- c("CallL2Fam<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
