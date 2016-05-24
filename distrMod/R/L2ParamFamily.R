## generating function
L2ParamFamily <- function(name, distribution = Norm(), distrSymm,
                          main = main(param), nuisance = nuisance(param),
                          fixed = fixed(param), trafo = trafo(param),
                          param = ParamFamParameter(name = paste("Parameter of", 
                                      name),  main = main, nuisance = nuisance, 
                                              fixed = fixed, trafo = trafo),
                          props = character(0),
                          startPar = NULL, makeOKPar = NULL,
                          modifyParam = function(theta){ Norm(mean=theta) },
                          L2deriv.fct = function(param) {force(theta <- param@main)
                                       return(function(x) {x-theta})},
                          L2derivSymm, L2derivDistr, L2derivDistrSymm,
                          FisherInfo.fct,
                          FisherInfo = FisherInfo.fct(param),
                          .returnClsName = NULL, .withMDE = TRUE){
     
    if(missing(name))
        name <- "L_2 differentiable parametric family of probability measures"
    if(missing(param)&&missing(main))
        param <- ParamFamParameter(name = "location", main = 0)
    if(missing(param)){
        argList <- list(name = paste("Parameter of", name),
                                   main = main)
        if(!missing(nuisance)) argList <- c(argList, nuisance = nuisance)                            
        if(!missing(fixed))    argList <- c(argList, fixed = fixed)                            
        if(!missing(trafo))    argList <- c(argList, trafo = trafo)                            
        param <- do.call(ParamFamParameter, argList)
        
    }
    if(missing(distrSymm)) distrSymm <- NoSymmetry()
    if(!is(distrSymm, "NoSymmetry")){
        if(!is(distrSymm@SymmCenter, "numeric"))
            stop("slot 'SymmCenter' of 'distrSymm' has to be of class 'numeric'")
        if(length(distrSymm@SymmCenter) != dimension(img(distribution)))
            stop("slot 'SymmCenter' of 'distrSymm' has wrong dimension")
    }
    fct <- L2deriv.fct(param)
    L2deriv0 <- if(!is.list(fct))
       RealRandVariable(list(fct), Domain = Reals()) else
       RealRandVariable(fct, Domain = Reals())
    L2deriv <- EuclRandVarList(L2deriv0)
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
    if(!is(L2derivDistr,"DistrList"))
         L2derivDistr <- UnivarDistrList(L2derivDistr)    
    if(!is(L2derivDistr,"UnivarDistrList"))
         L2derivDistr <- as(L2derivDistr,"UnivarDistrList")

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
                FI0 <- E(object = distribution, fun = L2deriv0 %*% t(L2deriv0 ))
        FisherInfo <- PosSemDefSymmMatrix(FI0)
    }else{
        FisherInfo <- PosSemDefSymmMatrix(FisherInfo)
    }
    if(ncol(FisherInfo) != dims)
        stop(paste("dimension of 'FisherInfo' should be", dims))

    if(missing(FisherInfo.fct))
        FisherInfo.fct <- function(param){
        fct <- L2deriv.fct(param)
        L2 <- if(!is.list(fct))
           RealRandVariable(list(fct), Domain = Reals()) else
           RealRandVariable(fct, Domain = Reals())
        return(PosSemDefSymmMatrix(E(object = distribution,
                                     fun = L2 %*% t(L2))))
        }

    parv <- c(param@main,param@nuisance)
    nms <- names(parv)
    
    if(!is.null(nms))
       dimnames(FisherInfo) <- list(nms,nms)

    f.call <- substitute(L2ParamFamily(name = N,
               distribution = D,
               distrSymm = DS,
               param = P,
               props = Props,
               startPar = sP,
               makeOKPar = okP,
               modifyParam = modP,
               L2deriv.fct = L2fct,
               L2derivSymm = L2Symm,
               L2derivDistr = L2D,
               L2derivDistrSymm = L2DSymm,
               FisherInfo.fct = Ffct,
               FisherInfo = FInfo,
               .returnClsName = rtn,
               .withMDE = wMDE0),
          list(N = name,
               D = distribution,
               DS = distrSymm,
               P = param,
               Props = props,
               sP = startPar,
               okP = makeOKPar,
               modP = modifyParam,
               L2fct = L2deriv.fct,
               L2Symm = L2derivSymm,
               L2D = L2derivDistr,
               L2DSymm = L2derivDistrSymm,
               Ffct = FisherInfo.fct,
               FInfo = FisherInfo,
               rtn = .returnClsName,
               wMDE0 =.withMDE))
 

    if(is.null(.returnClsName))
       .returnClsName <- "L2ParamFamily"
    L2Fam <- new(.returnClsName)
    L2Fam@name <- name
    L2Fam@distribution <- distribution
    L2Fam@fam.call <- f.call
    L2Fam@distrSymm <- distrSymm
    L2Fam@param <- param
    L2Fam@modifyParam <- modifyParam
    L2Fam@props <- props
    L2Fam@L2deriv.fct <- L2deriv.fct
    L2Fam@L2deriv <- L2deriv
    L2Fam@L2derivSymm <- L2derivSymm
    L2Fam@L2derivDistr <- L2derivDistr
    L2Fam@L2derivDistrSymm <- L2derivDistrSymm
    L2Fam@FisherInfo.fct <- FisherInfo.fct
    L2Fam@FisherInfo <- FisherInfo
    if(!is.null(startPar)) L2Fam@startPar <- startPar
    if(!is.null(makeOKPar)) L2Fam@makeOKPar <- makeOKPar
    L2Fam@.withMDE <- .withMDE
    return(L2Fam)
}

## access methods
setMethod("L2deriv", signature(object = "L2ParamFamily", param = "missing"), 
           function(object) object@L2deriv)
setMethod("L2deriv", signature(object = "L2ParamFamily", 
           param = "ParamFamParameter"), 
           function(object, param) object@L2deriv.fct(param))
setMethod("L2derivSymm", "L2ParamFamily", function(object) object@L2derivSymm)
setMethod("L2derivDistr", "L2ParamFamily", function(object){
                             ob <- object@L2derivDistr
                             if(is.call(ob)){
                                   ob <- eval(ob)
                                   eval.parent(object@L2derivDistr <- ob)
                             }
                             return(ob)
                             })
setMethod("L2derivDistrSymm", "L2ParamFamily", function(object) object@L2derivDistrSymm)
setMethod("FisherInfo", signature(object = "L2ParamFamily", param = "missing"),
           function(object) object@FisherInfo)
setMethod("FisherInfo", signature(object = "L2ParamFamily", param = "ParamFamParameter"),
           function(object, param) object@FisherInfo.fct(param))

## check centering of L2 derivative and Fisher Information
setMethod("checkL2deriv", "L2ParamFamily", 
    function(L2Fam, out = TRUE){ 
        dims <- length(L2Fam@param)
        L2deriv <- as(diag(dims) %*% L2Fam@L2deriv, "EuclRandVariable")

        cent <- E(object = L2Fam, fun = L2deriv)
        if(out) cat("precision of centering:\t", cent, "\n")

        consist <- E(object = L2Fam, fun = L2deriv %*% t(L2deriv))
        FI <- as(L2Fam@FisherInfo, "matrix")
        consist <- consist - FI
        if(out){
            cat("precision of Fisher information:\n")
            print(consist)
            cat("precision of Fisher information - relativ error [%]:\n")
            print(100*consist/FI)
        }

        if(out){
           cat("condition of Fisher information:\n")
           print(kappa(FI))
        }

        prec <- max(abs(cent), abs(consist))
 
        return(list(maximum.deviation = prec))
    })

