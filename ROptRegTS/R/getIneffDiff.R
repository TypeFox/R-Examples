###############################################################################
## get inefficiency difference for asymptotic MSE
###############################################################################
setMethod("getIneffDiff", signature(radius = "numeric", 
                                    L2Fam = "L2RegTypeFamily", 
                                    neighbor = "Neighborhood",
                                    risk = "asMSE"),
    function(radius, L2Fam, neighbor, risk, loRad, upRad, loRisk, upRisk, 
             z.start = NULL, A.start = NULL, upper.b, MaxIter, eps, warn){
        L2derivDim <- numberOfMaps(L2Fam@L2deriv)
        if(L2derivDim == 1){
            neighbor@radius <- radius
            res <- getInfRobRegTypeIC(ErrorL2deriv = L2Fam@ErrorL2derivDistr[[1]], 
                        Regressor = L2Fam@RegDistr, risk = risk, neighbor = neighbor, 
                        ErrorL2derivDistrSymm = L2Fam@ErrorL2derivDistrSymm[[1]], 
                        RegSymm = L2Fam@RegSymm, Finfo = L2Fam@FisherInfo, 
                        trafo = L2Fam@param@trafo, upper = upper.b, maxiter = MaxIter, 
                        tol = eps, warn = warn)
            trafo <- L2Fam@param@trafo
            ineffLo <- (sum(diag(res$A %*% t(trafo))) - res$b^2*(radius^2-loRad^2))/loRisk
            if(upRad == Inf)
                ineffUp <- res$b^2/upRisk
            else
                ineffUp <- (sum(diag(res$A %*% t(trafo))) - res$b^2*(radius^2-upRad^2))/upRisk
            assign("ineff", ineffUp, envir = sys.frame(which = -4))
            if(is(L2Fam@RegDistr, "MultivariateDistribution"))
                cat("current radius:\t", radius, "\tMSE-inefficiency difference:\t", ineffUp - ineffLo, "\n")

            return(ineffUp - ineffLo)
        }else{
            if(is(L2Fam@ErrorDistr, "UnivariateDistribution")){
                if((length(L2Fam@ErrorL2deriv) == 1) 
                   & is(L2Fam@ErrorL2deriv[[1]], "RealRandVariable")){
                    ErrorL2deriv <- L2Fam@ErrorL2deriv[[1]]
                    ErrorL2derivSymm <- L2Fam@ErrorL2derivSymm
                    ErrorL2derivDistrSymm <- L2Fam@ErrorL2derivDistrSymm
                }else{
                    ErrorL2deriv <- diag(dimension(L2Fam@ErrorL2deriv)) %*% L2Fam@ErrorL2deriv
                    ErrorL2deriv <- RealRandVariable(Map = ErrorL2deriv@Map, Domain = ErrorL2deriv@Domain)
                    nrvalues <- numberOfMaps(ErrorL2deriv)
                    if(numberOfMaps(L2Fam@ErrorL2deriv) != nrvalues){
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
                trafo <- L2Fam@param@trafo
                neighbor@radius <- radius
                res <- getInfRobRegTypeIC(ErrorL2deriv = ErrorL2deriv, 
                            Regressor = L2Fam@RegDistr, risk = risk, neighbor = neighbor, 
                            ErrorSymm = L2Fam@ErrorSymm, RegSymm = L2Fam@RegSymm, 
                            ErrorDistr = L2Fam@ErrorDistr, ErrorL2derivSymm = ErrorL2derivSymm, 
                            ErrorL2derivDistrSymm = ErrorL2derivDistrSymm,
                            Finfo = L2Fam@FisherInfo, trafo = trafo, upper = upper.b, 
                            z.start = z.start, A.start = A.start, maxiter = MaxIter, tol = eps, 
                            warn = warn)
                ineffLo <- (sum(diag(res$A%*%t(trafo))) - res$b^2*(radius^2-loRad^2))/loRisk
                if(upRad == Inf)
                    ineffUp <- res$b^2/upRisk
                else
                    ineffUp <- (sum(diag(res$A%*%t(trafo))) - res$b^2*(radius^2-upRad^2))/upRisk
                assign("ineff", ineffUp, envir = sys.frame(which = -4))
                cat("current radius:\t", radius, "\tMSE-inefficiency difference:\t", ineffUp - ineffLo, "\n")

                return(ineffUp - ineffLo)
            }else{
                stop("not yet implemented")
            }
        }
    })
setMethod("getIneffDiff", signature(radius = "numeric", 
                                    L2Fam = "L2RegTypeFamily", 
                                    neighbor = "Av2CondContNeighborhood",
                                    risk = "asMSE"),
    function(radius, L2Fam, neighbor, risk, loRad, upRad, loRisk, upRisk, 
             z.start = NULL, A.start = NULL, upper.b, MaxIter, eps, warn){
        L2derivDim <- numberOfMaps(L2Fam@L2deriv)
        if(L2derivDim == 1){
            neighbor@radius <- radius
            res <- getInfRobRegTypeIC(ErrorL2deriv = L2Fam@ErrorL2derivDistr[[1]], 
                        Regressor = L2Fam@RegDistr, risk = risk, neighbor = neighbor, 
                        ErrorL2derivDistrSymm = L2Fam@ErrorL2derivDistrSymm[[1]], 
                        RegSymm = L2Fam@RegSymm, Finfo = L2Fam@FisherInfo, 
                        trafo = L2Fam@param@trafo, upper = upper.b, maxiter = MaxIter, 
                        tol = eps, warn = warn)
            trafo <- L2Fam@param@trafo
            K.inv <- solve(E(L2Fam@RegDistr, fun = function(x){ x %*% t(x) }))
            ineffLo <- (res$A*sum(diag(t(trafo) %*% K.inv)) - res$b^2*(radius^2-loRad^2))/loRisk
            if(upRad == Inf)
                ineffUp <- res$b^2/upRisk
            else
                ineffUp <- (res$A*sum(diag(t(trafo) %*% K.inv)) - res$b^2*(radius^2-upRad^2))/upRisk
            assign("ineff", ineffUp, envir = sys.frame(which = -4))
#            cat("current radius:\t", radius, "\tMSE-inefficiency difference:\t", ineffUp - ineffLo, "\n")

            return(ineffUp - ineffLo)
        }else{
            stop("not yet implemented")
        }
    })
