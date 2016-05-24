###############################################################################
## radius minimax optimally robust IC 
## for L2ParamFamily and asymptotic risks
###############################################################################
setMethod("leastFavorableRadius", signature(L2Fam = "L2RegTypeFamily", 
                                            neighbor = "Neighborhood",
                                            risk = "asGRisk"),
    function(L2Fam, neighbor, risk, rho, upRad = 1, z.start = NULL, 
            A.start = NULL, upper = 100, maxiter = 100, 
            tol = .Machine$double.eps^0.4, warn = FALSE){
        ow <- options("warn")
        on.exit(options(ow))
        
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
                on.exit(options(ow))
                options(warn = -1)
                if(identical(all.equal(loRad, 0), TRUE)){
                    loRad <- 0
                    loRisk <- 1/as.vector(L2Fam@FisherInfo)
                }else{
                    neighbor@radius <- loRad
                    resLo <- getInfRobRegTypeIC(ErrorL2deriv = L2Fam@ErrorL2derivDistr[[1]], 
                                Regressor = L2Fam@RegDistr, risk = risk, neighbor = neighbor, 
                                ErrorL2derivDistrSymm = L2Fam@ErrorL2derivDistrSymm[[1]], 
                                RegSymm = L2Fam@RegSymm, Finfo = L2Fam@FisherInfo, 
                                trafo = L2Fam@param@trafo, upper = upper.b, maxiter = MaxIter, 
                                tol = eps, warn = warn)
                    loRisk <- getAsRiskRegTS(risk = risk, ErrorL2deriv = L2Fam@ErrorL2derivDistr[[1]], 
                                        Regressor = L2Fam@RegDistr, neighbor = neighbor, 
                                        clip = resLo$b, cent = resLo$a, stand = resLo$A, 
                                        trafo = L2Fam@param@trafo)[[1]]
                }

                if(upRad == Inf){
                    bmin <- getAsRiskRegTS(risk = asBias(), ErrorL2deriv = L2Fam@ErrorL2derivDistr[[1]], 
                                Regressor = L2Fam@Regressor, neighbor = neighbor, 
                                ErrorL2derivDistrSymm = L2Fam@ErrorL2derivDistrSymm[[1]], 
                                trafo = L2Fam@param@trafo, maxiter = MaxIter, tol = eps)$asBias
                    upRisk <- bmin^2
                }else{
                    neighbor@radius <- upRad
                    resUp <- getInfRobRegTypeIC(ErrorL2deriv = L2Fam@ErrorL2derivDistr[[1]], 
                                Regressor = L2Fam@RegDistr, risk = risk, neighbor = neighbor, 
                                ErrorL2derivDistrSymm = L2Fam@ErrorL2derivDistrSymm[[1]], 
                                RegSymm = L2Fam@RegSymm, Finfo = L2Fam@FisherInfo, 
                                trafo = L2Fam@param@trafo, upper = upper.b, maxiter = MaxIter, 
                                tol = eps, warn = warn)
                    upRisk <- getAsRiskRegTS(risk = risk, ErrorL2deriv = L2Fam@ErrorL2derivDistr[[1]], 
                                        Regressor = L2Fam@RegDistr, neighbor = neighbor, 
                                        clip = resUp$b, cent = resUp$a, stand = resUp$A, 
                                        trafo = L2Fam@param@trafo)[[1]]
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
                leastFavFct <- function(r, L2Fam, neighbor, risk, rho, 
                                        z.start, A.start, upper.b, MaxIter, eps, warn){
                    loRad <- r*rho
                    upRad <- r/rho
                    lower <- ifelse(identical(all.equal(loRad, 0), TRUE), 1e-4, loRad)
                    upper <- ifelse(upRad == Inf, 10, upRad)
                    ow <- options("warn")
                    on.exit(options(ow))
                    options(warn = -1)
                    trafo <- L2Fam@param@trafo
                    if(identical(all.equal(loRad, 0), TRUE)){
                        loRad <- 0
                        loRisk <- sum(diag(solve(L2Fam@FisherInfo)))
                    }else{
                        neighbor@radius <- loRad
                        resLo <- getInfRobRegTypeIC(ErrorL2deriv = ErrorL2deriv, 
                                    Regressor = L2Fam@RegDistr, risk = risk, neighbor = neighbor, 
                                    ErrorSymm = L2Fam@ErrorSymm, RegSymm = L2Fam@RegSymm, 
                                    ErrorDistr = L2Fam@ErrorDistr, ErrorL2derivSymm = ErrorL2derivSymm, 
                                    ErrorL2derivDistrSymm = ErrorL2derivDistrSymm,
                                    Finfo = L2Fam@FisherInfo, trafo = trafo, upper = upper.b, 
                                    z.start = z.start, A.start = A.start, maxiter = MaxIter, tol = eps, 
                                    warn = warn)
                        loRisk <- getAsRiskRegTS(risk = risk, ErrorL2deriv = ErrorL2deriv, 
                                            Regressor = L2Fam@RegDistr, neighbor = neighbor, 
                                            clip = resLo$b, cent = resLo$a, stand = resLo$A, 
                                            trafo = trafo)[[1]]
                    }

                    if(upRad == Inf){
                        bmin <- getAsRiskRegTS(risk = asBias(), ErrorL2deriv = ErrorL2deriv[[1]], 
                                    Regressor = L2Fam@Regressor, neighbor = neighbor, 
                                    ErrorDistr = L2Fam@ErrorDistr, trafo = trafo, 
                                    z.start = z.start, A.start = A.start, maxiter = MaxIter, 
                                    tol = eps)$asBias
                        upRisk <- bmin^2
                    }else{
                        neighbor@radius <- upRad
                        resUp <- getInfRobRegTypeIC(ErrorL2deriv = ErrorL2deriv, 
                                    Regressor = L2Fam@RegDistr, risk = risk, neighbor = neighbor, 
                                    ErrorSymm = L2Fam@ErrorSymm, RegSymm = L2Fam@RegSymm, 
                                    ErrorDistr = L2Fam@ErrorDistr, ErrorL2derivSymm = ErrorL2derivSymm, 
                                    ErrorL2derivDistrSymm = ErrorL2derivDistrSymm,
                                    Finfo = L2Fam@FisherInfo, trafo = trafo, upper = upper.b, 
                                    z.start = z.start, A.start = A.start, maxiter = MaxIter, tol = eps, 
                                    warn = warn)
                        upRisk <- getAsRiskRegTS(risk = risk, ErrorL2deriv = ErrorL2deriv, 
                                            Regressor = L2Fam@RegDistr, neighbor = neighbor, 
                                            clip = resUp$b, cent = resUp$a, stand = resUp$A, 
                                            trafo = trafo)[[1]]
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
                
                if(is.null(z.start)){ 
                    if(is(neighbor, "Av1CondContNeighborhood")){
                        k <- dimension(L2Fam@ErrorL2deriv)
                        z.start <- function(x){ numeric(k) }
                    }else{
                        z.start <- numeric(ncol(L2Fam@param@trafo))
                    }
                }

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
