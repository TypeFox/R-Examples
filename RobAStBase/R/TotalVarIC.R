## Generating function
TotalVarIC <- function(name, CallL2Fam = call("L2ParamFamily"),
                    Curve = EuclRandVarList(RealRandVariable(Map = c(function(x){x}), Domain = Reals())), 
                    Risks, Infos, clipLo = -Inf, clipUp = Inf, stand = as.matrix(1),
                    lowerCase = NULL, neighborRadius = 0, w = new("BdStWeight"),
                    normtype = NormType(), biastype = symmetricBias(),
                    modifyIC = NULL){

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
    if(!identical(dim(trafo(L2Fam@param)), dim(stand)))
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
    IC1@weight <- w
    IC1@biastype <- biastype
    IC1@normtype <- normtype
    IC1@modifyIC <- modifyIC

    return(IC1)
}

## generate IC
## for internal use only!
setMethod("generateIC", signature(neighbor = "TotalVarNeighborhood", 
                                  L2Fam = "L2ParamFamily"),
    function(neighbor, L2Fam, res){
        res$A <- A <- stand(res$w)
        clipLo <- clip(res$w)[1]
        clipUp <- clip(res$w)[2]
#        clipLo <- sign(as.vector(A))*res$a
        b <- clipUp-clipLo
        w <- res$w 
        if((clipLo == -Inf) & (b == Inf))
            clipUp <- Inf
        else
            clipUp <- clipLo + b

        L2call <- L2Fam@fam.call
        L2call$trafo <- trafo(L2Fam)

        cuv <- generateIC.fct(neighbor, L2Fam, res)
        return(TotalVarIC(
                name = "IC of total variation type", 
                CallL2Fam = L2call,
                Curve = generateIC.fct(neighbor, L2Fam, res),
                clipUp = clipUp,
                clipLo = clipLo,
                stand = A,
                lowerCase = res$d, 
                w = w,
                modifyIC = res$modifyIC,
                normtype = res$normtype,
                biastype = res$biastype,
                neighborRadius = neighbor@radius,
                Risks = res$risk,
                Infos = matrix(res$info, ncol = 2, 
                            dimnames = list(character(0), c("method", "message")))))
    })

## Access methods
setMethod("clip", "TotalVarIC", function(x1) x1@clipUp-x1@clipLo)
setMethod("clipLo", "TotalVarIC", function(object) object@clipLo)
setMethod("clipUp", "TotalVarIC", function(object) object@clipUp)
setMethod("neighbor", "TotalVarIC", function(object) TotalVarNeighborhood(radius = object@neighborRadius) )

## Replace methods
setReplaceMethod("clipLo", "TotalVarIC", 
    function(object, value){ 
        stopifnot(is.numeric(value))
        L2Fam <- eval(object@CallL2Fam)
        w <- object@weight
        clip(w)[1] <- value
        weight(w) <- getweight(w, neighbor = TotalVarNeighborhood(radius = object@neighborRadius), 
                               biastype = object@biastype, 
                               normW = object@normtype)
        res <- list(A = object@stand, a = value, b = object@clipUp-value, 
                    d = object@lowerCase, risk = object@Risks, info = object@Infos, 
                    w = w, normtype = object@normtype, biastype = object@biastype,
                    modifyIC = object@modifyIC)
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
        w <- object@weight
        clip(w)[2] <- value
        weight(w) <- getweight(w, neighbor = TotalVarNeighborhood(radius = object@neighborRadius), 
                               biastype = object@biastype, 
                               normW = object@normtype)
        res <- list(A = object@stand, a = object@clipLo, b = value-object@clipLo, 
                    d = object@lowerCase, risk = object@Risks, info = object@Infos, 
                    w = w, normtype = object@normtype, biastype = object@biastype,
                    modifyIC = object@modifyIC)
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
        w <- object@weight
        stand(w) <- value
        weight(w) <- getweight(w, neighbor = TotalVarNeighborhood(radius = object@neighborRadius), 
                               biastype = object@biastype, 
                               normW = object@normtype)
        res <- list(A = value, a = object@clipLo, b = object@clipUp-object@clipLo, 
                    d = object@lowerCase, risk = object@Risks, info = object@Infos, 
                    w = w, normtype = object@normtype, biastype = object@biastype,
                    modifyIC = object@modifyIC)
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
                    d = value, risk = object@Risks, info = object@Infos, 
                    w = object@weight, normtype = object@normtype, biastype = object@biastype,
                    modifyIC = object@modifyIC)
        object <- generateIC(neighbor = TotalVarNeighborhood(radius = object@neighborRadius), 
                             L2Fam = L2Fam, res = res)
        addInfo(object) <- c("lowerCase<-", "The slot 'lowerCase' has been changed")
        addInfo(object) <- c("lowerCase<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
setReplaceMethod("CallL2Fam", "TotalVarIC",
    function(object, value){ 
        L2Fam <- eval(value)
        res <- list(A = object@stand, a = object@clipLo, b = object@clipUp-object@clipLo, 
                    d = object@lowerCase, risk = object@Risks, info = object@Infos, 
                    w = object@weight, normtype = object@normtype, biastype = object@biastype,
                    modifyIC = object@modifyIC)
        object <- generateIC(neighbor = TotalVarNeighborhood(radius = object@neighborRadius), 
                             L2Fam = L2Fam, res = res)
        addInfo(object) <- c("CallL2Fam<-", "The slot 'CallL2Fam' has been changed")
        addInfo(object) <- c("CallL2Fam<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
