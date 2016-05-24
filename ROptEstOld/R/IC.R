## generating function
IC <- function(name, Curve = EuclRandVarList(RealRandVariable(Map = list(function(x){x}), 
               Domain = Reals())), Risks, Infos, CallL2Fam = call("L2ParamFamily")){
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
    trafo <- L2Fam@param@trafo
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
    
    return(IC1)
}

## access methods
setMethod("CallL2Fam", "IC", function(object) object@CallL2Fam)

## replace methods
setReplaceMethod("CallL2Fam", "IC",
    function(object, value){
        object@CallL2Fam <- value
        object
    })

## check centering and Fisher consistency
setMethod("checkIC", signature(IC = "IC", L2Fam = "missing"), 
    function(IC, out = TRUE){ 
        L2Fam <- eval(IC@CallL2Fam)
        trafo <- L2Fam@param@trafo
        IC1 <- as(diag(nrow(trafo)) %*% IC@Curve, "EuclRandVariable")
        
        cent <- E(L2Fam, IC1)
        if(out)
            cat("precision of centering:\t", cent, "\n")

        dims <- length(L2Fam@param)
        L2deriv <- as(diag(dims) %*% L2Fam@L2deriv, "EuclRandVariable")
        consist <- E(L2Fam, IC1 %*% t(L2deriv)) - trafo
        if(out){
            cat("precision of Fisher consistency:\n")
            print(consist)
        }
        prec <- max(abs(cent), abs(consist))
        names(prec) <- "maximum deviation"
        
        return(prec)
    })
## check centering and Fisher consistency
setMethod("checkIC", signature(IC = "IC", L2Fam = "L2ParamFamily"), 
    function(IC, L2Fam, out = TRUE){ 
        D1 <- L2Fam@distribution
        if(dimension(Domain(IC@Curve[[1]])) != dimension(img(D1)))
            stop("dimension of 'Domain' of 'Curve' != dimension of 'img' of 'distribution' of 'L2Fam'")

        trafo <- L2Fam@param@trafo
        IC1 <- as(diag(dimension(IC@Curve)) %*% IC@Curve, "EuclRandVariable")
        cent <- E(D1, IC1)
        if(out)
            cat("precision of centering:\t", cent, "\n")

        dims <- length(L2Fam@param)
        L2deriv <- as(diag(dims) %*% L2Fam@L2deriv, "EuclRandVariable")

        consist <- E(D1, IC1 %*% t(L2deriv)) - trafo
        if(out){
            cat("precision of Fisher consistency:\n")
            print(consist)
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
