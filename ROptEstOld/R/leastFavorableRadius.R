###############################################################################
## radius minimax optimally robust IC 
## for L2ParamFamily and asymptotic risks
###############################################################################
setMethod("leastFavorableRadius", signature(L2Fam = "L2ParamFamily", 
                                            neighbor = "UncondNeighborhood",
                                            risk = "asGRisk"),
    function(L2Fam, neighbor, risk, rho, upRad = 1, z.start = NULL, 
            A.start = NULL, upper = 100, maxiter = 100, 
            tol = .Machine$double.eps^0.4, warn = FALSE){
        if(length(rho) != 1)
            stop("'rho' is not of length == 1")
        if((rho <= 0)||(rho >= 1))
            stop("'rho' not in (0,1)")

        L2derivDim <- numberOfMaps(L2Fam@L2deriv)
        if(L2derivDim == 1){
            leastFavFct <- function(r, L2Fam, neighbor, risk, rho, 
                                    upper.b, MaxIter, eps, warn){
                loRad <- r*rho
                upRad <- r/rho
                lower <- ifelse(identical(all.equal(loRad, 0), TRUE), 1e-4, loRad)
                upper <- ifelse(upRad == Inf, 10, upRad)
                ow <- options("warn")
                options(warn = -1)
                if(identical(all.equal(loRad, 0), TRUE)){
                    loRad <- 0
                    loRisk <- 1/as.vector(L2Fam@FisherInfo)
                }else{
                    neighbor@radius <- loRad
                    resLo <- getInfRobIC(L2deriv = L2Fam@L2derivDistr[[1]], neighbor = neighbor, 
                                risk = risk, symm = L2Fam@L2derivDistrSymm[[1]],
                                Finfo = L2Fam@FisherInfo, upper = upper.b,
                                trafo = L2Fam@param@trafo, maxiter = MaxIter, tol = eps, warn = warn)
                    loRisk <- getAsRisk(risk = risk, L2deriv = L2Fam@L2derivDistr[[1]], 
                                        neighbor = neighbor, clip = resLo$b, cent = resLo$a, 
                                        stand = resLo$A, trafo = L2Fam@param@trafo)[[1]]
                }

                if(upRad == Inf){
                    bmin <- getAsRisk(risk = asBias(), L2deriv = L2Fam@L2derivDistr[[1]], 
                                neighbor = neighbor, trafo = L2Fam@param@trafo, symm = L2Fam@L2derivSymm[[1]])
                    upRisk <- bmin^2
                }else{
                    neighbor@radius <- upRad
                    resUp <- getInfRobIC(L2deriv = L2Fam@L2derivDistr[[1]], neighbor = neighbor, 
                                risk = risk, symm = L2Fam@L2derivDistrSymm[[1]],
                                Finfo = L2Fam@FisherInfo, upper = upper.b,
                                trafo = L2Fam@param@trafo, maxiter = MaxIter, tol = eps, warn = warn)
                    upRisk <- getAsRisk(risk = risk, L2deriv = L2Fam@L2derivDistr[[1]], 
                                        neighbor = neighbor, clip = resUp$b, cent = resUp$a, 
                                        stand = resUp$A, trafo = L2Fam@param@trafo)[[1]]
                }
                leastFavR <- uniroot(getIneffDiff, lower = lower, upper = upper, 
                                tol = .Machine$double.eps^0.25, L2Fam = L2Fam, neighbor = neighbor, 
                                risk = risk, loRad = loRad, upRad = upRad, loRisk = loRisk, 
                                upRisk = upRisk, upper.b = upper.b, eps = eps, MaxIter = MaxIter, 
                                warn = warn)$root
                options(ow)
                cat("current radius:\t", r, "\tinefficiency:\t", ineff, "\n")
                return(ineff)
            }
            leastFavR <- optimize(leastFavFct, lower = 1e-4, upper = upRad, 
                            tol = .Machine$double.eps^0.25, maximum = TRUE,
                            L2Fam = L2Fam, neighbor = neighbor, risk = risk,
                            rho = rho, upper.b = upper, MaxIter = maxiter, 
                            eps = tol, warn = warn)

            res <- list(rho = rho, leastFavorableRadius = leastFavR$maximum, 
                        ineff = leastFavR$objective)
            names(res)[3] <- paste(class(risk)[1], "-inefficiency", sep="")
            return(res)
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
                leastFavFct <- function(r, L2Fam, neighbor, risk, rho, 
                                        z.start, A.start, upper.b, MaxIter, eps, warn){
                    loRad <- r*rho
                    upRad <- r/rho
                    lower <- ifelse(identical(all.equal(loRad, 0), TRUE), 1e-4, loRad)
                    upper <- ifelse(upRad == Inf, 10, upRad)
                    ow <- options("warn")
                    options(warn = -1)
                    trafo <- L2Fam@param@trafo
                    if(identical(all.equal(loRad, 0), TRUE)){
                        loRad <- 0
                        loRisk <- sum(diag(solve(L2Fam@FisherInfo)))
                    }else{
                        neighbor@radius <- loRad
                        resLo <- getInfRobIC(L2deriv = L2deriv, neighbor = neighbor, risk = risk, 
                                    Distr = L2Fam@distribution, DistrSymm = L2Fam@distrSymm, 
                                    L2derivSymm = L2derivSymm, L2derivDistrSymm = L2derivDistrSymm, 
                                    Finfo = L2Fam@FisherInfo, trafo = trafo, z.start = z.start, 
                                    A.start = A.start, upper = upper.b, maxiter = MaxIter, 
                                    tol = eps, warn = warn)
                        loRisk <- getAsRisk(risk = risk, L2deriv = L2deriv, neighbor = neighbor, 
                                        clip = resLo$b, cent = resLo$a, stand = resLo$A, trafo = trafo)[[1]]
                    }

                    if(upRad == Inf){
                        bmin <- getAsRisk(risk = asBias(), L2deriv = L2deriv, neighbor = neighbor, 
                                    Distr = L2Fam@distribution, L2derivDistrSymm = L2Fam@L2derivDistrSymm,
                                    trafo = trafo, z.start = z.start, A.start = A.start, 
                                    maxiter = maxiter, tol = tol)$asBias
                        upRisk <- bmin^2
                    }else{
                        neighbor@radius <- upRad
                        resUp <- getInfRobIC(L2deriv = L2deriv, neighbor = neighbor, risk = risk, 
                                    Distr = L2Fam@distribution, DistrSymm = L2Fam@distrSymm, 
                                    L2derivSymm = L2derivSymm, L2derivDistrSymm = L2derivDistrSymm, 
                                    Finfo = L2Fam@FisherInfo, trafo = trafo, z.start = z.start, 
                                    A.start = A.start, upper = upper.b, maxiter = maxiter, 
                                    tol = tol, warn = warn)
                         upRisk <- getAsRisk(risk = risk, L2deriv = L2deriv, neighbor = neighbor, 
                                        clip = resUp$b, cent = resUp$a, stand = resUp$A, trafo = trafo)[[1]]
                    }
                    leastFavR <- uniroot(getIneffDiff, lower = lower, upper = upper, 
                                    tol = .Machine$double.eps^0.25, L2Fam = L2Fam, neighbor = neighbor, 
                                    z.start = z.start, A.start = A.start, upper.b = upper.b, risk = risk, 
                                    loRad = loRad, upRad = upRad, loRisk = loRisk, upRisk = upRisk,
                                    eps = eps, MaxIter = MaxIter, warn = warn)$root
                    options(ow)
                    cat("current radius:\t", r, "\tinefficiency:\t", ineff, "\n")
                    return(ineff)
                }
                if(is.null(z.start)) z.start <- numeric(L2derivDim)
                if(is.null(A.start)) A.start <- L2Fam@param@trafo
                leastFavR <- optimize(leastFavFct, lower = 1e-4, upper = upRad, 
                                tol = .Machine$double.eps^0.25, maximum = TRUE,
                                L2Fam = L2Fam, neighbor = neighbor, risk = risk,
                                rho = rho, z.start = z.start, A.start = A.start, 
                                upper.b = upper, MaxIter = maxiter, eps = tol, warn = warn)

                res <- list(rho = rho, leastFavorableRadius = leastFavR$maximum, 
                            ineff = leastFavR$objective)
                names(res)[3] <- paste(class(risk)[1], "-inefficiency", sep="")
                return(res)
            }else{
                stop("not yet implemented")
            }
        }
    })
