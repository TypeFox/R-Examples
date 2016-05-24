###############################################################################
## radius minimax optimally robust IC 
## for L2ParamFamily and asymptotic risks
###############################################################################
setMethod("radiusMinimaxIC", signature(L2Fam = "L2ParamFamily", 
                                       neighbor = "UncondNeighborhood",
                                       risk = "asGRisk"),
    function(L2Fam, neighbor, risk, loRad, upRad, z.start = NULL, A.start = NULL,
            upper = 1e5, maxiter = 100, tol = .Machine$double.eps^0.4, warn = FALSE){
        if(length(loRad) != 1)
            stop("'loRad' is not of length == 1")
        if(length(upRad) != 1)
            stop("'upRad' is not of length == 1")
        if(loRad >= upRad)
            stop("'upRad < loRad' is not fulfilled")
        L2derivDim <- numberOfMaps(L2Fam@L2deriv)
        if(L2derivDim == 1){
            ow <- options("warn")
            options(warn = -1)
            upper.b <- upper
            lower <- ifelse(identical(all.equal(loRad, 0), TRUE), 1e-4, loRad)
            upper <- ifelse(upRad == Inf, max(loRad+1, 2), upRad)

            if(identical(all.equal(loRad, 0), TRUE)){
                loRad <- 0
                loRisk <- 1/as.vector(L2Fam@FisherInfo)
            }else{
                neighbor@radius <- loRad
                resLo <- getInfRobIC(L2deriv = L2Fam@L2derivDistr[[1]], neighbor = neighbor, 
                            risk = risk, symm = L2Fam@L2derivDistrSymm[[1]],
                            Finfo = L2Fam@FisherInfo, upper = upper.b,
                            trafo = L2Fam@param@trafo, maxiter = maxiter, tol = tol, warn = warn)
                loRisk <- getAsRisk(risk = risk, L2deriv = L2Fam@L2derivDistr[[1]], 
                                    neighbor = neighbor, clip = resLo$b, cent = resLo$a, 
                                    stand = resLo$A, trafo = L2Fam@param@trafo)[[1]]
            }

            if(upRad == Inf){
                bmin <- getAsRisk(risk = asBias(), L2deriv = L2Fam@L2derivDistr[[1]], 
                            neighbor = neighbor, trafo = L2Fam@param@trafo)$asBias
                upRisk <- bmin^2
            }else{
                neighbor@radius <- upRad
                resUp <- getInfRobIC(L2deriv = L2Fam@L2derivDistr[[1]], neighbor = neighbor, 
                            risk = risk, symm = L2Fam@L2derivDistrSymm[[1]],
                            Finfo = L2Fam@FisherInfo, upper = upper.b,
                            trafo = L2Fam@param@trafo, maxiter = maxiter, tol = tol, warn = warn)
                upRisk <- getAsRisk(risk = risk, L2deriv = L2Fam@L2derivDistr[[1]], 
                                    neighbor = neighbor, clip = resUp$b, cent = resUp$a, 
                                    stand = resUp$A, trafo = L2Fam@param@trafo)[[1]]
            }

            leastFavR <- uniroot(getIneffDiff, lower = lower, upper = upper, 
                            tol = .Machine$double.eps^0.25, L2Fam = L2Fam, neighbor = neighbor, 
                            upper.b = upper.b, risk = risk, loRad = loRad, upRad = upRad, 
                            loRisk = loRisk, upRisk = upRisk, eps = tol, 
                            MaxIter = maxiter, warn = warn)$root
            neighbor@radius <- leastFavR
            res <- getInfRobIC(L2deriv = L2Fam@L2derivDistr[[1]], neighbor = neighbor, 
                        risk = risk, symm = L2Fam@L2derivSymm[[1]],
                        Finfo = L2Fam@FisherInfo, upper = upper.b,
                        trafo = L2Fam@param@trafo, maxiter = maxiter, tol = tol, warn = warn)
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
                ow <- options("warn")
                options(warn = -1)
                upper.b <- upper
                lower <- ifelse(identical(all.equal(loRad, 0), TRUE), 1e-4, loRad)
                upper <- ifelse(upRad == Inf, max(loRad+1, 2), upRad)

                if(identical(all.equal(loRad, 0), TRUE)){
                    loRad <- 0
                    loRisk <- sum(diag(solve(L2Fam@FisherInfo)))
                }else{
                    neighbor@radius <- loRad
                    resLo <- getInfRobIC(L2deriv = L2deriv, neighbor = neighbor, risk = risk, 
                                Distr = L2Fam@distribution, DistrSymm = L2Fam@distrSymm, 
                                L2derivSymm = L2derivSymm, L2derivDistrSymm = L2derivDistrSymm, 
                                Finfo = L2Fam@FisherInfo, trafo = trafo, z.start = z.start, 
                                A.start = A.start, upper = upper.b, maxiter = maxiter, 
                                tol = tol, warn = warn)
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
                                eps = tol, MaxIter = maxiter, warn = warn)$root
                neighbor@radius <- leastFavR
                res <- getInfRobIC(L2deriv = L2deriv, neighbor = neighbor, risk = risk, 
                            Distr = L2Fam@distribution, DistrSymm = L2Fam@distrSymm, 
                            L2derivSymm = L2derivSymm, L2derivDistrSymm = L2derivDistrSymm, 
                            Finfo = L2Fam@FisherInfo, trafo = trafo, z.start = z.start, 
                            A.start = A.start, upper = upper.b, maxiter = maxiter, 
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
                stop("not yet implemented")
            }
        }
    })
