## generating function
L2ParamFamily <- function(name, distribution = Norm(), distrSymm, 
                          main = 0, nuisance, trafo, param, props = character(0), 
                          L2deriv = EuclRandVarList(RealRandVariable(list(function(x){x}), Domain = Reals())),
                          L2derivSymm, L2derivDistr, L2derivDistrSymm, FisherInfo){
    if(missing(name)) 
        name <- "L_2 differentiable parametric family of probability measures"
    if(missing(distrSymm)) distrSymm <- NoSymmetry()
    if(!is(distrSymm, "NoSymmetry")){
        if(!is(distrSymm@SymmCenter, "numeric"))
            stop("slot 'SymmCenter' of 'distrSymm' has to be of class 'numeric'")
        if(length(distrSymm@SymmCenter) != dimension(img(distribution)))
            stop("slot 'SymmCenter' of 'distrSymm' has wrong dimension")
    }
    if(missing(param)) 
        param <- ParamFamParameter(name = paste("Parameter of", name), 
                                   main = main, nuisance = nuisance, trafo = trafo)
    if(missing(L2derivSymm)){
        nrvalues <- numberOfMaps(L2deriv)
        L <- vector("list", nrvalues)
        for(i in 1:nrvalues) L[[i]] <- NonSymmetric()
        L2derivSymm <- new("FunSymmList", L)
    }
    if(is(distribution, "UnivariateCondDistribution"))
        stop("conditional distributions are not allowed in slot 'distribution'")
    if(missing(L2derivDistr))
        L2derivDistr <- imageDistr(RandVar = L2deriv, distr = distribution)
    if(missing(L2derivDistrSymm)){
        nrvalues <- length(L2derivDistr)
        L <- vector("list", nrvalues)
        for(i in 1:nrvalues) L[[i]] <- NoSymmetry()
        L2derivDistrSymm <- new("DistrSymmList", L)
    }

    nrvalues <- numberOfMaps(L2deriv)
    if(nrvalues != length(L2derivSymm))
        stop("number of Maps of 'L2deriv' != length of 'L2derivSymm'")
    if(nrvalues != length(L2derivDistr))
        stop("number of Maps of 'L2deriv' != length of 'L2derivDistr'")
    if(nrvalues != length(L2derivDistrSymm))
        stop("number of Maps of 'L2deriv' != length of 'L2derivDistrSymm'")
    if(dimension(Domain(L2deriv[[1]])) != dimension(img(distribution)))
        stop("dimension of 'Domain' of 'L2deriv' != dimension of 'img' of 'distribution'")
    dims <- length(param)
    if(dimension(L2deriv) != dims)
        stop("dimension of 'L2deriv' != dimension of parameters")

    if(missing(FisherInfo)){
        L2 <- as(diag(dims) %*% L2deriv, "EuclRandVariable")        
        FisherInfo <- PosDefSymmMatrix(E(object = distribution, fun = L2 %*% t(L2)))
    }else{
        FisherInfo <- PosDefSymmMatrix(FisherInfo)
    }
    if(ncol(FisherInfo) != dims)
        stop(paste("dimension of 'FisherInfo' should be", dims))
    
    L2Fam <- new("L2ParamFamily")
    L2Fam@name <- name
    L2Fam@distribution <- distribution
    L2Fam@distrSymm <- distrSymm
    L2Fam@param <- param
    L2Fam@props <- props
    L2Fam@L2deriv <- L2deriv
    L2Fam@L2derivSymm <- L2derivSymm
    L2Fam@L2derivDistr <- L2derivDistr
    L2Fam@L2derivDistrSymm <- L2derivDistrSymm
    L2Fam@FisherInfo <- FisherInfo

    return(L2Fam)
}

## access methods
setMethod("L2deriv", "L2ParamFamily", function(object) object@L2deriv)
setMethod("L2derivSymm", "L2ParamFamily", function(object) object@L2derivSymm)
setMethod("L2derivDistr", "L2ParamFamily", function(object) object@L2derivDistr)
setMethod("L2derivDistrSymm", "L2ParamFamily", function(object) object@L2derivDistrSymm)
setMethod("FisherInfo", "L2ParamFamily", function(object) object@FisherInfo)

## check centering of L2 derivative and Fisher Information
setMethod("checkL2deriv", "L2ParamFamily", 
    function(L2Fam, out = TRUE){ 
        dims <- length(L2Fam@param)
        L2deriv <- as(diag(dims) %*% L2Fam@L2deriv, "EuclRandVariable")

        cent <- E(object = L2Fam, fun = L2deriv)
        if(out) cat("precision of centering:\t", cent, "\n")
        
        consist <- E(object = L2Fam, fun = L2deriv %*% t(L2deriv))
        consist <- consist - as(L2Fam@FisherInfo, "matrix")
        if(out){
            cat("precision of Fisher information:\n")
            print(consist)
        }
        
        prec <- max(abs(cent), abs(consist))
        
        return(list(maximum.deviation = prec))
    })
