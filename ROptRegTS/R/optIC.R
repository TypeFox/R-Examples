###############################################################################
## Classical optimal IC (optimal in sense of the Cramer-Rao bound)
###############################################################################
setMethod("optIC", signature(model = "L2RegTypeFamily", risk = "asCov"),
    function(model, risk){
        Curve <- as((model@param@trafo %*% solve(model@FisherInfo)) %*% model@L2deriv, "EuclRandVariable")
        return(IC(
            name = paste("Classical optimal influence curve for", model@name), 
            CallL2Fam = call("L2RegTypeFamily", 
                            name = model@name,
                            distribution = model@distribution,
                            distrSymm = model@distrSymm,  
                            param = model@param,
                            props = model@props,
                            L2deriv = model@L2deriv,
                            ErrorDistr = model@ErrorDistr,
                            ErrorSymm = model@ErrorSymm,
                            RegDistr = model@RegDistr,
                            RegSymm = model@RegSymm,
                            Regressor = model@Regressor,
                            ErrorL2deriv = model@ErrorL2deriv,
                            ErrorL2derivSymm = model@ErrorL2derivSymm,
                            ErrorL2derivDistr = model@ErrorL2derivDistr,
                            ErrorL2derivDistrSymm = model@ErrorL2derivDistrSymm,
                            FisherInfo = model@FisherInfo),
            Curve = EuclRandVarList(Curve), 
            Risks = list(asCov = model@param@trafo %*% solve(model@FisherInfo) %*% t(model@param@trafo)),
            Infos = matrix(c("optIC", "optimal IC in sense of Cramer-Rao bound"), 
                        ncol = 2, dimnames = list(character(0), c("method", "message")))))
    })

###############################################################################
## Optimally robust IC for infinitesimal robust regression type model 
## and asymptotic risks
###############################################################################
setMethod("optIC", signature(model = "InfRobRegTypeModel", risk = "asRisk"),
    function(model, risk, z.start=NULL, A.start=NULL, upper = 1e4, 
             maxiter = 50, tol = .Machine$double.eps^0.4, warn = TRUE){
        ErrorL2derivDim <- numberOfMaps(model@center@ErrorL2deriv)
        ow <- options("warn")
        on.exit(options(ow))
        if(ErrorL2derivDim == 1){
            options(warn = -1)
            res <- getInfRobRegTypeIC(ErrorL2deriv = model@center@ErrorL2derivDistr[[1]], 
                        Regressor = model@center@RegDistr, risk = risk, neighbor = model@neighbor, 
                        ErrorL2derivDistrSymm = model@center@ErrorL2derivDistrSymm[[1]], 
                        RegSymm = model@center@RegSymm, Finfo = model@center@FisherInfo, 
                        trafo = model@center@param@trafo, upper = upper, 
                        maxiter = maxiter, tol = tol, warn = warn)
            options(ow)             
            res$info <- c("optIC", res$info)
            return(generateIC(model@neighbor, model@center, res))
        }else{
            if(is(model@center@ErrorDistr, "UnivariateDistribution")){
                if((length(model@center@ErrorL2deriv) == 1) 
                   & is(model@center@ErrorL2deriv[[1]], "RealRandVariable")){
                    ErrorL2deriv <- model@center@ErrorL2deriv[[1]]
                    ErrorL2derivSymm <- model@center@ErrorL2derivSymm
                    ErrorL2derivDistrSymm <- model@center@ErrorL2derivDistrSymm
                }else{
                    ErrorL2deriv <- diag(dimension(model@center@ErrorL2deriv)) %*% model@center@ErrorL2deriv
                    ErrorL2deriv <- RealRandVariable(Map = ErrorL2deriv@Map, Domain = ErrorL2deriv@Domain)
                    nrvalues <- numberOfMaps(ErrorL2deriv)
                    if(numberOfMaps(model@center@ErrorL2deriv) != nrvalues){
                        L1 <- vector("list", nrvalues)
                        L2 <- vector("list", nrvalues)
                        for(i in 1:nrvalues){
                            L1[[i]] <- NonSymmetric()
                            L2[[i]] <- NoSymmetry()
                        }
                        ErrorL2derivSymm <- new("FunSymmList", L1)
                        ErrorL2derivDistrSymm <- new("DistrSymmList", L2)
                    }
                }
                options(warn = -1)
                res <- getInfRobRegTypeIC(ErrorL2deriv = ErrorL2deriv, 
                            Regressor = model@center@RegDistr, risk = risk, neighbor = model@neighbor, 
                            ErrorSymm = model@center@ErrorSymm, 
                            RegSymm = model@center@RegSymm, ErrorDistr = model@center@ErrorDistr,
                            ErrorL2derivSymm = ErrorL2derivSymm, ErrorL2derivDistrSymm = ErrorL2derivDistrSymm,
                            Finfo = model@center@FisherInfo, trafo = model@center@param@trafo, 
                            upper = upper, z.start = z.start, A.start = A.start, maxiter = maxiter, 
                            tol = tol, warn = warn)
                options(ow)                   
                res$info <- c("optIC", res$info)
                return(generateIC(model@neighbor, model@center, res))
            }else{
                stop("not yet implemented")
            }
        }
    })
###############################################################################
## Optimally robust IC for infinitesimal robust regression type model 
## and asymptotic under-/overshoot risk
###############################################################################
setMethod("optIC", signature(model = "InfRobRegTypeModel", risk = "asUnOvShoot"),
    function(model, risk, upper = 1e4, maxiter = 50, tol = .Machine$double.eps^0.4, 
             warn = TRUE){
        ow <- options("warn")
        on.exit(options(ow))
        ErrorL2derivDistr <- model@center@ErrorL2derivDistr[[1]]
        if((length(model@center@ErrorL2derivDistr) == 1) & is(ErrorL2derivDistr, "UnivariateDistribution")
            & is(model@center@RegDistr, "UnivariateDistribution")){
            options(warn = -1)
            res <- getInfRobRegTypeIC(ErrorL2deriv = ErrorL2derivDistr, 
                        Regressor = model@center@RegDistr, risk = risk, neighbor = model@neighbor, 
                        ErrorL2derivDistrSymm = model@center@ErrorL2derivDistrSymm[[1]], 
                        RegSymm = model@center@RegSymm, Finfo = model@center@FisherInfo, 
                        trafo = model@center@param@trafo, upper = upper, 
                        maxiter = maxiter, tol = tol, warn = warn)
            options(ow)
            if(is(model@neighbor, "UncondNeighborhood")){
                if(is(model@neighbor, "ContNeighborhood"))
                    res$info <- c("optIC", "optIC", res$info, "Optimal IC for 'InfRobRegTypeModel' with 'ContNeighborhood'!!!")
                else
                    res$info <- c("optIC", res$info)
                return(generateIC(TotalVarNeighborhood(radius = model@neighbor@radius), model@center, res))
            }else{
                if(is(model@neighbor, "CondContNeighborhood"))
                    res$info <- c("optIC", "optIC", res$info, "Optimal IC for 'InfRobRegTypeModel' with 'CondContNeighborhood'!!!")
                else
                    res$info <- c("optIC", res$info)
                return(generateIC(CondTotalVarNeighborhood(radius = model@neighbor@radius, 
                                                           radiusCurve = model@neighbor@radiusCurve), 
                                  model@center, res))
            }
        }else{
            stop("restricted to linear regression with 1-dimensional regressors")
        }
    })
###############################################################################
## Optimally robust IC for fixed robust regression type model 
## and finite-sample under-/overshoot risk
###############################################################################
setMethod("optIC", signature(model = "FixRobRegTypeModel", risk = "fiUnOvShoot"),
    function(model, risk, sampleSize, upper = 1e4, maxiter = 50, 
             tol = .Machine$double.eps^0.4, warn = TRUE, Algo = "A", cont = "left"){
        ow <- options("warn")
        on.exit(options(ow))
        if(!identical(all.equal(sampleSize, trunc(sampleSize)), TRUE))
            stop("'sampleSize' has to be an integer > 0")
        if(is(model@center@ErrorDistr, "UnivariateDistribution") 
           & is(model@center@RegDistr, "UnivariateDistribution")){
            RegDistr <- model@center@RegDistr
            if(!is(RegDistr, "AbscontDistribution"))
                if(!identical(all.equal(d(RegDistr)(0), 0), TRUE))
                    stop("Solution only available under 'K(x=0)!=0'!") 
            options(warn = -1)
            res <- getFixRobRegTypeIC(ErrorDistr = model@center@ErrorDistr, 
                        Regressor = RegDistr, risk = risk, neighbor = model@neighbor, 
                        sampleSize = sampleSize, upper = upper, maxiter = maxiter, 
                        tol = tol, warn = warn, Algo = Algo, cont = cont)
            options(ow)
            if(is(model@neighbor, "UncondNeighborhood")){
                if(is(model@neighbor, "ContNeighborhood"))
                    res$info <- c("optIC", "optIC", res$info, "Optimal IC for 'FixRobRegTypeModel' with 'ContNeighborhood'!!!")
                else
                    res$info <- c("optIC", res$info)
                return(generateIC(TotalVarNeighborhood(radius = model@neighbor@radius), model@center, res))
            }else{
                if(is(model@neighbor, "CondContNeighborhood"))
                    res$info <- c("optIC", "optIC", res$info, "Optimal IC for 'FixRobRegTypeModel' with 'CondContNeighborhood'!!!")
                else
                    res$info <- c("optIC", res$info)
                return(generateIC(CondTotalVarNeighborhood(radius = model@neighbor@radius, 
                                                           radiusCurve = model@neighbor@radiusCurve), 
                                  model@center, res))
            }
        }else{
            stop("restricted to linear regression with 1-dimensional regressors")
        }
    })
