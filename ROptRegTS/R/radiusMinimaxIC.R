###############################################################################
## radius minimax optimally robust IC 
## for L2ParamFamily and asymptotic risks
###############################################################################
setMethod("radiusMinimaxIC", signature(L2Fam = "L2RegTypeFamily", 
                                       neighbor = "Neighborhood",
                                       risk = "asGRisk"),
    function(L2Fam, neighbor, risk, loRad, upRad, z.start = NULL, A.start = NULL,
            upper = 1e4, maxiter = 100, tol = .Machine$double.eps^0.4, warn = FALSE){
        
        ow <- options("warn")
        on.exit(options(ow))
        
        if(length(loRad) != 1)
            stop("'loRad' is not of length == 1")
        if(length(upRad) != 1)
            stop("'upRad' is not of length == 1")
        if(loRad >= upRad)
            stop("'upRad < loRad' is not fulfilled")
        L2derivDim <- numberOfMaps(L2Fam@ErrorL2deriv)
        if(L2derivDim == 1){
            options(warn = -1)
            upper.b <- upper
            lower <- ifelse(identical(all.equal(loRad, 0), TRUE), 1e-4, loRad)
            upper <- ifelse(upRad == Inf, max(loRad+1, 2), upRad)

            if(identical(all.equal(loRad, 0), TRUE)){
                loRad <- 0
                loRisk <- sum(diag(solve(L2Fam@FisherInfo)))
            }else{
                neighbor@radius <- loRad
                resLo <- getInfRobRegTypeIC(ErrorL2deriv = L2Fam@ErrorL2derivDistr[[1]], 
                            Regressor = L2Fam@RegDistr, risk = risk, neighbor = neighbor, 
                            ErrorL2derivDistrSymm = L2Fam@ErrorL2derivDistrSymm[[1]], 
                            RegSymm = L2Fam@RegSymm, Finfo = L2Fam@FisherInfo, 
                            trafo = L2Fam@param@trafo, upper = upper.b, maxiter = maxiter, 
                            tol = tol, warn = warn)
                loRisk <- getAsRiskRegTS(risk = risk, ErrorL2deriv = L2Fam@ErrorL2derivDistr[[1]], 
                                    Regressor = L2Fam@RegDistr, neighbor = neighbor, 
                                    clip = resLo$b, cent = resLo$a, stand = resLo$A, 
                                    trafo = L2Fam@param@trafo)[[1]]
            }

            if(upRad == Inf){
                bmin <- getAsRiskRegTS(risk = asBias(), ErrorL2deriv = L2Fam@ErrorL2derivDistr[[1]], 
                            Regressor = L2Fam@RegDistr, neighbor = neighbor, 
                            ErrorL2derivDistrSymm = L2Fam@ErrorL2derivDistrSymm[[1]], 
                            trafo = L2Fam@param@trafo, maxiter = maxiter, tol = tol)$asBias
                upRisk <- bmin^2
            }else{
                neighbor@radius <- upRad
                resUp <- getInfRobRegTypeIC(ErrorL2deriv = L2Fam@ErrorL2derivDistr[[1]], 
                            Regressor = L2Fam@RegDistr, risk = risk, neighbor = neighbor, 
                            ErrorL2derivDistrSymm = L2Fam@ErrorL2derivDistrSymm[[1]], 
                            RegSymm = L2Fam@RegSymm, Finfo = L2Fam@FisherInfo, 
                            trafo = L2Fam@param@trafo, upper = upper.b, maxiter = maxiter, 
                            tol = tol, warn = warn)
                upRisk <- getAsRiskRegTS(risk = risk, ErrorL2deriv = L2Fam@ErrorL2derivDistr[[1]], 
                                    Regressor = L2Fam@RegDistr, neighbor = neighbor, 
                                    clip = resUp$b, cent = resUp$a, stand = resUp$A, 
                                    trafo = L2Fam@param@trafo)[[1]]
            }

            leastFavR <- uniroot(getIneffDiff, lower = lower, upper = upper, 
                            tol = .Machine$double.eps^0.25, L2Fam = L2Fam, neighbor = neighbor, 
                            upper.b = upper.b, risk = risk, loRad = loRad, upRad = upRad, 
                            loRisk = loRisk, upRisk = upRisk, eps = tol, 
                            MaxIter = maxiter, warn = warn)$root
            neighbor@radius <- leastFavR
            res <- getInfRobRegTypeIC(ErrorL2deriv = L2Fam@ErrorL2derivDistr[[1]], 
                        Regressor = L2Fam@RegDistr, risk = risk, neighbor = neighbor, 
                        ErrorL2derivDistrSymm = L2Fam@ErrorL2derivDistrSymm[[1]], 
                        RegSymm = L2Fam@RegSymm, Finfo = L2Fam@FisherInfo, 
                        trafo = L2Fam@param@trafo, upper = upper.b, maxiter = maxiter, 
                        tol = tol, warn = warn)
            options(ow)                   
            res$info <- c("radiusMinimaxIC", paste("radius minimax IC for radius interval [", 
                            round(loRad, 3), ", ", round(upRad, 3), "]", sep=""))
            res$info <- rbind(res$info, c("radiusMinimaxIC", 
                            paste("least favorable radius: ", round(leastFavR, 3), sep="")))
            res$info <- rbind(res$info, c("radiusMinimaxIC", 
                            paste("maximum ", sQuote(class(risk)[1]), "-inefficiency: ",
                            round(ineff, 3), sep="")))
            return(generateIC(neighbor, L2Fam, res))
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
                options(warn = -1)
                upper.b <- upper
                lower <- ifelse(identical(all.equal(loRad, 0), TRUE), 1e-4, loRad)
                upper <- ifelse(upRad == Inf, max(loRad+1, 2), upRad)

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
                                z.start = z.start, A.start = A.start, maxiter = maxiter, tol = tol, 
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
                                z.start = z.start, A.start = A.start, maxiter = maxiter, 
                                tol = tol)$asBias
                    upRisk <- bmin^2
                }else{
                    neighbor@radius <- upRad
                    resUp <- getInfRobRegTypeIC(ErrorL2deriv = ErrorL2deriv, 
                                Regressor = L2Fam@RegDistr, risk = risk, neighbor = neighbor, 
                                ErrorSymm = L2Fam@ErrorSymm, RegSymm = L2Fam@RegSymm, 
                                ErrorDistr = L2Fam@ErrorDistr, ErrorL2derivSymm = ErrorL2derivSymm, 
                                ErrorL2derivDistrSymm = ErrorL2derivDistrSymm,
                                Finfo = L2Fam@FisherInfo, trafo = trafo, upper = upper.b, 
                                z.start = z.start, A.start = A.start, maxiter = maxiter, tol = tol, 
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
                                eps = tol, MaxIter = maxiter, warn = warn)$root
                neighbor@radius <- leastFavR
                res <- getInfRobRegTypeIC(ErrorL2deriv = ErrorL2deriv, 
                            Regressor = L2Fam@RegDistr, risk = risk, neighbor = neighbor, 
                            ErrorSymm = L2Fam@ErrorSymm, RegSymm = L2Fam@RegSymm, 
                            ErrorDistr = L2Fam@ErrorDistr, ErrorL2derivSymm = ErrorL2derivSymm, 
                            ErrorL2derivDistrSymm = ErrorL2derivDistrSymm,
                            Finfo = L2Fam@FisherInfo, trafo = trafo, upper = upper.b, 
                            z.start = z.start, A.start = A.start, maxiter = maxiter, tol = tol, 
                            warn = warn)
                options(ow)                   
                res$info <- c("radiusMinimaxIC", paste("radius minimax IC for radius interval [", 
                                round(loRad, 3), ", ", round(upRad, 3), "]", sep=""))
                res$info <- rbind(res$info, c("radiusMinimaxIC", 
                                paste("least favorable radius: ", round(leastFavR, 3), sep="")))
                res$info <- rbind(res$info, c("radiusMinimaxIC", 
                                paste("maximum ", sQuote(class(risk)[1]), "-inefficiency: ",
                            round(ineff, 3), sep="")))
                return(generateIC(neighbor, L2Fam, res))
            }else{
                stop("not yet implemented")
            }
        }
    })
