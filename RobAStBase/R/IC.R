## generating function
IC <- function(name, Curve = EuclRandVarList(RealRandVariable(Map = list(function(x){x}), 
               Domain = Reals())), Risks, Infos, CallL2Fam = call("L2ParamFamily"),
               modifyIC = NULL){
    if(missing(name))
        name <- "square integrable (partial) influence curve"
    if(missing(Risks))
        Risks <- list()
    if(missing(Infos))
        Infos <- matrix(c(character(0),character(0)), ncol=2,
                     dimnames=list(character(0), c("method", "message")))

    if(!is(Domain(Curve[[1]]), "EuclideanSpace"))
        stop("The domain of 'Curve' has to be a Euclidean space")
    if(!is.character(Infos))
        stop("'Infos' contains no matrix of characters")
    for(char in names(Risks))
        if(!extends(char, "RiskType"))
            stop(paste(char, "is no valid 'RiskType'"))
    if(ncol(Infos)!=2)
        stop("'Infos' must have two columns")

    L2Fam <- eval(CallL2Fam)
    trafo <- trafo(L2Fam@param)
    if(nrow(trafo) != dimension(Curve))
        stop("wrong dimension of 'Curve'")
    if(dimension(Domain(L2Fam@L2deriv[[1]])) != dimension(Domain(Curve[[1]])))
        stop("dimension of 'Domain' of 'L2deriv' != dimension of 'Domain' of 'Curve'")

    IC1 <- new("IC")
    IC1@name <- name
    IC1@Curve <- Curve
    IC1@Risks <- Risks
    IC1@Infos <- Infos
    IC1@CallL2Fam <- CallL2Fam
    IC1@modifyIC <- modifyIC

    return(IC1)
}

## access methods
setMethod("CallL2Fam", "IC", function(object) object@CallL2Fam)
setMethod("modifyIC", "IC", function(object) object@modifyIC)

## replace methods
setReplaceMethod("CallL2Fam", "IC",
    function(object, value){
        object@CallL2Fam <- value
        object
    })

## check centering and Fisher consistency
setMethod("checkIC", signature(IC = "IC", L2Fam = "missing"), 
    function(IC, out = TRUE, ...){
        L2Fam <- eval(IC@CallL2Fam)
        checkIC(IC, L2Fam, out = out, ...)
    })
## check centering and Fisher consistency
setMethod("checkIC", signature(IC = "IC", L2Fam = "L2ParamFamily"), 
    function(IC, L2Fam, out = TRUE, ...){
        D1 <- L2Fam@distribution
        if(dimension(Domain(IC@Curve[[1]])) != dimension(img(D1)))
            stop("dimension of 'Domain' of 'Curve' != dimension of 'img' of 'distribution' of 'L2Fam'")

        trafo <- trafo(L2Fam@param)
        IC1 <- as(diag(dimension(IC@Curve)) %*% IC@Curve, "EuclRandVariable")
        cent <- E(D1, IC1, ...)
        if(out)
            cat("precision of centering:\t", cent, "\n")

        dims <- length(L2Fam@param)
        L2deriv <- as(diag(dims) %*% L2Fam@L2deriv, "EuclRandVariable")

        consist <- E(D1, IC1 %*% t(L2deriv), ...) - trafo
        if(out){
            cat("precision of Fisher consistency:\n")
            print(consist)
            cat("precision of Fisher consistency - relativ error [%]:\n")
            print(100*consist/trafo)
        }

        prec <- max(abs(cent), abs(consist))
        names(prec) <- "maximum deviation"

        return(prec)
    })

## evaluate IC
setMethod("evalIC", signature(IC = "IC", x = "numeric"), 
    function(IC, x){ 
        if(!is.null(IC@Curve[[1]]@Domain)){
            if(length(x) != IC@Curve[[1]]@Domain@dimension)
                stop("x has wrong dimension")
        }

        dimn <- dimension(IC@Curve)
        Curve <- as(diag(dimn) %*% IC@Curve, "EuclRandVariable")

        return(as.vector(evalRandVar(Curve, x)))
    })
setMethod("evalIC", signature(IC = "IC", x = "matrix"), 
    function(IC, x){ 
        if(!is.null(IC@Curve[[1]]@Domain)){
            if(ncol(x) != IC@Curve[[1]]@Domain@dimension)
                stop("x has wrong dimension")
        }

        dimn <- dimension(IC@Curve)
        Curve <- as(diag(dimn) %*% IC@Curve, "EuclRandVariable")

        if(dimn == 1)
            return(t(evalRandVar(Curve, x)[,,1]))
        else
            return(evalRandVar(Curve, x)[,,1])
    })

## make some L2function a pIC at a model
setMethod("makeIC", signature(IC = "IC", L2Fam = "missing"), 
    function(IC){ 
        L2Fam <- eval(IC@CallL2Fam)
        makeIC(IC, L2Fam)
    })

## make some L2function a pIC at a model
setMethod("makeIC", signature(IC = "IC", L2Fam = "L2ParamFamily"), 
    function(IC, L2Fam){ 
        D1 <- L2Fam@distribution
        if(dimension(Domain(IC@Curve[[1]])) != dimension(img(D1)))
            stop("dimension of 'Domain' of 'Curve' != dimension of 'img' of 'distribution' of 'L2Fam'")

        trafo <- trafo(L2Fam@param)
        IC1 <- as(diag(dimension(IC@Curve)) %*% IC@Curve, "EuclRandVariable")
        cent <- E(D1, IC1)
        IC1 <- IC1 - cent

        dims <- length(L2Fam@param)
        if(dimension(IC@Curve) != dims)
           stop("Dimension of IC and parameter must be equal")

        L2deriv <- as(diag(dims) %*% L2Fam@L2deriv, "EuclRandVariable")

        E1 <- matrix(E(L2Fam, IC1 %*% t(L2deriv)), dims, dims)
        stand <- trafo %*% solve(E1) 
        Y <- as(stand %*% IC1, "EuclRandVariable")
        #ICfct <- vector(mode = "list", length = dims)
        #ICfct[[1]] <- function(x){Y(x)}

        modifyIC <- function(L2Fam, IC){ makeIC(IC, L2Fam) }

        CallL2Fam <- L2Fam@fam.call

        return(IC(name = name(IC),
                  Curve = EuclRandVarList(Y),
                  Risks = list(), 
                  Infos=matrix(c("IC<-", 
                                 "generated by affine linear trafo to enforce consistency"), 
                               ncol=2, dimnames=list(character(0), c("method", "message"))), 
                  CallL2Fam = CallL2Fam,
                  modifyIC = modifyIC))
    })


# alias to IC needed here:
.IC <- IC

setMethod("makeIC", signature(IC = "list", L2Fam = "L2ParamFamily"),
    function(IC, L2Fam, forceIC = TRUE, name, Risks, Infos, modifyIC = NULL){
        mc <- match.call(call = sys.call(sys.parent(1)), expand.dots = FALSE)[-1]
        mc0 <- as.list(mc)
        mc0$IC <- NULL
        mc0$L2Fam <- NULL
        mc0$forceIC <- NULL
        if(!all(as.logical(c(lapply(IC,is.function)))))
           stop("First argument must be a list of functions")

        IC.1 <- lapply(IC, function(IC.2) 
                  if(length(formals(IC.2))==0) function(x) IC.2(x) else IC.2)

        mc0$Curve <- EuclRandVarList(RealRandVariable(Map = IC.1, Domain = Reals()))
        mc0$CallL2Fam <- substitute(L2Fam@fam.call)

        IC.0 <- do.call(.IC,mc0)
        if(forceIC) IC.0 <- makeIC(IC.0, L2Fam)
        return(IC.0)
    })



setMethod("makeIC", signature(IC = "function", L2Fam = "L2ParamFamily"),
    function(IC, L2Fam, forceIC = TRUE, name, Risks, Infos, modifyIC = NULL){
        mc <- match.call(call = sys.call(sys.parent(1)), expand.dots = FALSE)[-1]
        mc0 <- as.list(mc)
        mc0$IC <- NULL
        mc0$L2Fam <- NULL
        mc0$forceIC <- NULL
        IC.1 <- if(length(formals(IC))==0) function(x) IC(x) else IC
        mc0$Curve <- EuclRandVarList(RealRandVariable(Map = list(IC.1),
                         Domain = Reals()))
        mc0$CallL2Fam <- substitute(L2Fam@fam.call)
        print(mc0)
        
        IC.0 <- do.call(.IC,mc0)
        print(IC.0)
        if(forceIC) IC.0 <- makeIC(IC.0, L2Fam)
        return(IC.0)
    })
