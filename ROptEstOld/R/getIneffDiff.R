###############################################################################
## get inefficiency difference for asymptotic MSE
###############################################################################
setMethod("getIneffDiff", signature(radius = "numeric", 
                                    L2Fam = "L2ParamFamily", 
                                    neighbor = "UncondNeighborhood",
                                    risk = "asMSE"),
    function(radius, L2Fam, neighbor, risk, loRad, upRad, loRisk, upRisk, 
             z.start = NULL, A.start = NULL, upper.b, MaxIter, eps, warn){
        L2derivDim <- numberOfMaps(L2Fam@L2deriv)
        if(L2derivDim == 1){
            neighbor@radius <- radius
            res <- getInfRobIC(L2deriv = L2Fam@L2derivDistr[[1]], neighbor = neighbor, 
                        risk = risk, symm = L2Fam@L2derivDistrSymm[[1]], 
                        Finfo = L2Fam@FisherInfo, upper = upper.b,
                        trafo = L2Fam@param@trafo, maxiter = MaxIter, tol = eps, warn = warn)
            trafo <- as.vector(L2Fam@param@trafo)
            ineffLo <- (as.vector(res$A)*trafo - res$b^2*(radius^2-loRad^2))/loRisk
            if(upRad == Inf)
                ineffUp <- res$b^2/upRisk
            else
                ineffUp <- (as.vector(res$A)*trafo - res$b^2*(radius^2-upRad^2))/upRisk
            assign("ineff", ineffUp, envir = sys.frame(which = -4))

            return(ineffUp - ineffLo)
        }else{
            if(is(L2Fam@distribution, "UnivariateDistribution")){
                if((length(L2Fam@L2deriv) == 1) & is(L2Fam@L2deriv[[1]], "RealRandVariable")){
                    L2deriv <- L2Fam@L2deriv[[1]]
                    L2derivSymm <- L2Fam@L2derivSymm
                    L2derivDistrSymm <- L2Fam@L2derivDistrSymm
                }else{
                    L2deriv <- diag(dimension(L2Fam@L2deriv)) %*% L2Fam@L2deriv
                    L2deriv <- RealRandVariable(Map = L2deriv@Map, Domain = L2deriv@Domain)
                    nrvalues <- numberOfMaps(L2deriv)
                    if(numberOfMaps(L2Fam@L2deriv) != nrvalues){
                        L1 <- vector("list", nrvalues)
                        L2 <- vector("list", nrvalues)
                        for(i in 1:nrvalues){
                            L1[[i]] <- NonSymmetric()
                            L2[[i]] <- NoSymmetry()
                        }
                        L2derivSymm <- new("FunSymmList", L1)
                        L2derivDistrSymm <- new("DistrSymmList", L2)
                    }
                }
                trafo <- L2Fam@param@trafo
                neighbor@radius <- radius
                res <- getInfRobIC(L2deriv = L2deriv, neighbor = neighbor, risk = risk, 
                            Distr = L2Fam@distribution, DistrSymm = L2Fam@distrSymm, 
                            L2derivSymm = L2derivSymm, L2derivDistrSymm = L2derivDistrSymm, 
                            Finfo = L2Fam@FisherInfo, trafo = trafo, z.start = z.start, 
                            A.start = A.start, upper = upper.b, maxiter = MaxIter, 
                            tol = eps, warn = warn)
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
