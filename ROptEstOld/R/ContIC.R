## Generating function
ContIC <- function(name, CallL2Fam = call("L2ParamFamily"),
                   Curve = EuclRandVarList(RealRandVariable(Map = c(function(x){x}), Domain = Reals())), 
                   Risks, Infos, clip = Inf, cent = 0, stand = as.matrix(1), 
                   lowerCase = NULL, neighborRadius = 0){
    if(missing(name))
        name <- "IC of contamination type"
    if(missing(Risks))
        Risks <- list()
    if(missing(Infos))
        Infos <- matrix(c(character(0),character(0)), ncol=2,
                    dimnames=list(character(0), c("method", "message")))

    if(any(neighborRadius < 0)) # radius vector?!
        stop("'neighborRadius' has to be in [0, Inf]")
    if(length(cent) != nrow(stand))
        stop("length of centering constant != nrow of standardizing matrix")
    if((length(clip) != 1) && (length(clip) != length(Curve)))
        stop("length of clipping bound != 1 and != length of 'Curve'")
    if(!is.null(lowerCase))
        if(length(lowerCase) != nrow(stand))
            stop("length of 'lowerCase' != nrow of standardizing matrix")
    L2Fam <- eval(CallL2Fam)
    if(!identical(dim(L2Fam@param@trafo), dim(stand)))
        stop(paste("dimension of 'trafo' of 'param' != dimension of 'stand'"))
 
    contIC <- new("ContIC")
    contIC@name <- name
    contIC@Curve <- Curve
    contIC@Risks <- Risks
    contIC@Infos <- Infos
    contIC@CallL2Fam <- CallL2Fam
    contIC@clip <- clip
    contIC@cent <- cent
    contIC@stand <- stand
    contIC@lowerCase <- lowerCase
    contIC@neighborRadius <- neighborRadius

    return(contIC)   
#    return(new("ContIC", name = name, Curve = Curve, Risks = Risks, Infos = Infos,
#               CallL2Fam = CallL2Fam, clip = clip, cent = cent, stand = stand, 
#               lowerCase = lowerCase, neighborRadius = neighborRadius))
}

## generate IC
## for internal use only!
setMethod("generateIC", signature(neighbor = "ContNeighborhood", 
                                  L2Fam = "L2ParamFamily"),
    function(neighbor, L2Fam, res){
        A <- res$A
        a <- res$a
        b <- res$b
        d <- res$d
        nrvalues <- nrow(A)
        ICfct <- vector(mode = "list", length = nrvalues)
        Y <- as(A %*% L2Fam@L2deriv - a, "EuclRandVariable")
        if(nrvalues == 1){
            if(!is.null(d)){
                ICfct[[1]] <- function(x){ 
                                    ind <- (Y(x) != 0) 
                                    b*(ind*Y(x)/(ind*absY(x) + (1-ind)) + zi*(1-ind)*d)
                              }
                body(ICfct[[1]]) <- substitute(
                                        { ind <- (Y(x) != 0) 
                                          b*(ind*Y(x)/(ind*absY(x) + (1-ind)) + zi*(1-ind)*d) },
                                        list(Y = Y@Map[[1]], absY = abs(Y)@Map[[1]], b = b, d = d, 
                                             zi = sign(L2Fam@param@trafo)))
            }else{
                ICfct[[1]] <- function(x){ Y(x)*pmin(1, b/absY(x)) }
                body(ICfct[[1]]) <- substitute({ Y(x)*pmin(1, b/absY(x)) },
                                                 list(Y = Y@Map[[1]], absY = abs(Y)@Map[[1]], b = b))
            }
        }else{
            absY <- sqrt(Y %*% Y)
            if(!is.null(d))
                for(i in 1:nrvalues){
                    ICfct[[i]] <- function(x){ ind <- (Yi(x) != 0) ; ind*b*Yi(x)/absY(x) + (1-ind)*d }
                    body(ICfct[[i]]) <- substitute({ ind <- (Yi(x) != 0) ; ind*b*Yi(x)/absY(x) + (1-ind)*d },
                                                 list(Yi = Y@Map[[i]], absY = absY@Map[[1]], b = b, d = d[i]))
                }
            else
                for(i in 1:nrvalues){
                    ICfct[[i]] <- function(x){ Yi(x)*pmin(1, b/absY(x)) }
                    body(ICfct[[i]]) <- substitute({ Yi(x)*pmin(1, b/absY(x)) },
                                                 list(Yi = Y@Map[[i]], absY = absY@Map[[1]], b = b))
                }
        }
        return(ContIC(
                name = "IC of contamination type", 
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
                clip = b,
                cent = a,
                stand = A,
                lowerCase = d,
                neighborRadius = neighbor@radius,
                Risks = res$risk,
                Infos = matrix(res$info, ncol = 2, 
                            dimnames = list(character(0), c("method", "message")))))
    })

## Access methods
setMethod("clip", "ContIC", function(object) object@clip)
setMethod("cent", "ContIC", function(object) object@cent)
setMethod("stand", "ContIC", function(object) object@stand)
setMethod("lowerCase", "ContIC", function(object) object@lowerCase)
setMethod("neighborRadius", "ContIC", function(object) object@neighborRadius)

## replace methods
setReplaceMethod("clip", "ContIC", 
    function(object, value){ 
        stopifnot(is.numeric(value))
        L2Fam <- eval(object@CallL2Fam)
        res <- list(A = object@stand, a = object@cent, b = value, d = object@lowerCase,
                    risk = object@Risks, info = object@Infos)
        object <- generateIC(neighbor = ContNeighborhood(radius = object@neighborRadius), 
                             L2Fam = L2Fam, res = res)
        addInfo(object) <- c("clip<-", "The clipping bound has been changed")
        addInfo(object) <- c("clip<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
setReplaceMethod("cent", "ContIC", 
    function(object, value){ 
        stopifnot(is.numeric(value))
        L2Fam <- eval(object@CallL2Fam)
        res <- list(A = object@stand, a = value, b = object@clip, d = object@lowerCase,
                    risk = object@Risks, info = object@Infos)
        object <- generateIC(neighbor = ContNeighborhood(radius = object@neighborRadius), 
                             L2Fam = L2Fam, res = res)
        addInfo(object) <- c("cent<-", "The centering constant has been changed")
        addInfo(object) <- c("cent<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
setReplaceMethod("stand", "ContIC", 
    function(object, value){ 
        stopifnot(is.matrix(value))
        L2Fam <- eval(object@CallL2Fam)
        res <- list(A = value, a = object@cent, b = object@clip, d = object@lowerCase,
                    risk = object@Risks, info = object@Infos)
        object <- generateIC(neighbor = ContNeighborhood(radius = object@neighborRadius), 
                             L2Fam = L2Fam, res = res)
        addInfo(object) <- c("stand<-", "The standardizing matrix has been changed")
        addInfo(object) <- c("stand<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
setReplaceMethod("lowerCase", "ContIC", 
    function(object, value){ 
        stopifnot(is.null(value)||is.numeric(value))
        L2Fam <- eval(object@CallL2Fam)
        res <- list(A = object@stand, a = object@cent, b = object@clip, d = value,
                    risk = object@Risks, info = object@Infos)
        object <- generateIC(neighbor = ContNeighborhood(radius = object@neighborRadius), 
                             L2Fam = L2Fam, res = res)
        addInfo(object) <- c("lowerCase<-", "The slot 'lowerCase' has been changed")
        addInfo(object) <- c("lowerCase<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
setReplaceMethod("neighborRadius", "ContIC", 
    function(object, value){ 
        object@neighborRadius <- value
        if(any(value < 0)) # radius vector?!
            stop("'value' has to be in [0, Inf]")
        addInfo(object) <- c("neighborRadius<-", "The slot 'neighborRadius' has been changed")
        addInfo(object) <- c("neighborRadius<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
setReplaceMethod("CallL2Fam", "ContIC",
    function(object, value){ 
        L2Fam <- eval(value)
        res <- list(A = object@stand, a = object@cent, b = object@clip, d = object@lowerCase,
                    risk = object@Risks, info = object@Infos)
        object <- generateIC(neighbor = ContNeighborhood(radius = object@neighborRadius), 
                             L2Fam = L2Fam, res = res)
        addInfo(object) <- c("CallL2Fam<-", "The slot 'CallL2Fam' has been changed")
        addInfo(object) <- c("CallL2Fam<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
