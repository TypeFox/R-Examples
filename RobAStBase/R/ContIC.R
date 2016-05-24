## Generating function
ContIC <- function(name, CallL2Fam = call("L2ParamFamily"),
                   Curve = EuclRandVarList(RealRandVariable(Map = c(function(x){x}), 
                                           Domain = Reals())), 
                   Risks, Infos, clip = Inf, cent = 0, stand = as.matrix(1), 
                   lowerCase = NULL, neighborRadius = 0, w = new("HampelWeight"),
                   normtype = NormType(), biastype = symmetricBias(),
                   modifyIC = NULL){
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
    if(!identical(dim(trafo(L2Fam@param)), dim(stand)))
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
    contIC@weight <- w
    contIC@biastype <- biastype
    contIC@normtype <- normtype
    contIC@modifyIC <- modifyIC

    return(contIC)
#    return(new("ContIC", name = name, Curve = Curve, Risks = Risks, Infos = Infos,
#               CallL2Fam = CallL2Fam, clip = clip, cent = cent, stand = stand, 
#               lowerCase = lowerCase, neighborRadius = neighborRadius))
}


setMethod("generateIC", signature(neighbor = "ContNeighborhood", 
                                  L2Fam = "L2ParamFamily"),
    function(neighbor, L2Fam, res){
        A <- res$A
        a <- res$a
        b <- res$b
        d <- res$d
        normtype <- res$normtype
        biastype <- res$biastype
        w <- res$w
        L2call <- L2Fam@fam.call
        L2call$trafo <- trafo(L2Fam)
        return(ContIC(
                name = "IC of contamination type", 
                CallL2Fam = L2call,
                Curve = generateIC.fct(neighbor, L2Fam, res),
                clip = b,
                cent = a,
                stand = A,
                lowerCase = d,
                w = w,
                neighborRadius = neighbor@radius,
                modifyIC = res$modifyIC,
                normtype = normtype,
                biastype = biastype,
                Risks = res$risk,
                Infos = matrix(res$info, ncol = 2, 
                            dimnames = list(character(0), c("method", "message")))))
    })

## Access methods
setMethod("clip", "ContIC", function(x1) x1@clip)
setMethod("cent", "ContIC", function(object) object@cent)
setMethod("neighbor", "ContIC", function(object) ContNeighborhood(radius = object@neighborRadius) )

## replace methods
setReplaceMethod("clip", "ContIC", 
    function(object, value){ 
        stopifnot(is.numeric(value))
        L2Fam <- eval(object@CallL2Fam)
        w <- object@weight
        clip(w) <- value
        weight(w) <- getweight(w, neighbor = ContNeighborhood(radius = object@neighborRadius), 
                               biastype = object@biastype, 
                               normW = object@normtype)
        res <- list(A = object@stand, a = object@cent, b = value, d = object@lowerCase,
                    risk = object@Risks, info = object@Infos, w = w,
                    normtype = object@normtype, biastype = object@biastype,
                    modifyIC = object@modifyIC)
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
        w <- object@weight
        cent(w) <- as.vector(solve(object@stand) %*% value)
        weight(w) <- getweight(w, neighbor = ContNeighborhood(radius = object@neighborRadius), 
                               biastype = object@biastype, 
                               normW = object@normtype)
        res <- list(A = object@stand, a = value, b = object@clip, d = object@lowerCase,
                    risk = object@Risks, info = object@Infos, w = w,
                    normtype = object@normtype, biastype = object@biastype,
                    modifyIC = object@modifyIC)
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
        w <- object@weight
        stand(w) <- value
        weight(w) <- getweight(w, neighbor = ContNeighborhood(radius = object@neighborRadius), 
                               biastype = object@biastype, 
                               normW = object@normtype)
        res <- list(A = value, a = object@cent, b = object@clip, d = object@lowerCase,
                    risk = object@Risks, info = object@Infos, w = w,
                    normtype = object@normtype, biastype = object@biastype,
                    modifyIC = object@modifyIC)
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
                    risk = object@Risks, info = object@Infos, w = object@weight,
                    normtype = object@normtype, biastype = object@biastype,
                    modifyIC = object@modifyIC)
        object <- generateIC(neighbor = ContNeighborhood(radius = object@neighborRadius), 
                             L2Fam = L2Fam, res = res)
        addInfo(object) <- c("lowerCase<-", "The slot 'lowerCase' has been changed")
        addInfo(object) <- c("lowerCase<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
setReplaceMethod("CallL2Fam", "ContIC",
    function(object, value){ 
        L2Fam <- eval(value)
        res <- list(A = object@stand, a = object@cent, b = object@clip, d = object@lowerCase,
                    risk = object@Risks, info = object@Infos, w = object@weight,
                    normtype = object@normtype, biastype = object@biastype,
                    modifyIC = object@modifyIC)
        object <- generateIC(neighbor = ContNeighborhood(radius = object@neighborRadius), 
                             L2Fam = L2Fam, res = res)
        addInfo(object) <- c("CallL2Fam<-", "The slot 'CallL2Fam' has been changed")
        addInfo(object) <- c("CallL2Fam<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
