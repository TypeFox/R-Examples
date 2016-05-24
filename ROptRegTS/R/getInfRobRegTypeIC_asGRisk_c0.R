###############################################################################
## get optimally robust IC for convex asymptotic risks
###############################################################################
setMethod("getInfRobRegTypeIC", signature(ErrorL2deriv = "UnivariateDistribution",
                                          Regressor = "Distribution", 
                                          risk = "asGRisk", 
                                          neighbor = "ContNeighborhood"),
    function(ErrorL2deriv, Regressor, risk, neighbor, ErrorL2derivDistrSymm, 
             RegSymm, Finfo, trafo, upper, maxiter, tol, warn){
        radius <- neighbor@radius
        if(identical(all.equal(radius, 0), TRUE)){
            if(warn) cat("'radius == 0' => (classical) optimal IC\n", 
                         "in sense of Cramer-Rao bound is returned\n")
            res <- getInfRobRegTypeIC(ErrorL2deriv = ErrorL2deriv, Regressor = Regressor,
                        risk = asCov(), neighbor = neighbor, 
                        ErrorL2derivDistrSymm = ErrorL2derivDistrSymm,
                        RegSymm = RegSymm, Finfo = Finfo, trafo = trafo)
            Risk <- getAsRiskRegTS(risk = risk, ErrorL2deriv = ErrorL2deriv, 
                                   Regressor = Regressor, neighbor = neighbor,
                                   clip = res$b, cent = res$a, stand = res$A, 
                                   trafo = trafo)
            res$risk <- c(Risk, res$risk)
            return(res)
        }
        z <- numeric(ncol(trafo))
        A <- trafo
        b <- 0
        if(is(ErrorL2derivDistrSymm, "SphericalSymmetry")) 
            z.comp <- !(ErrorL2derivDistrSymm@SymmCenter == 0)
        else
            z.comp <- TRUE

        iter <- 0
        repeat{
            iter <- iter + 1
            b.old <- b
            z.old <- z
            A.old <- A
            b <- try(uniroot(getInfClipRegTS, lower = .Machine$double.eps^0.75, 
                        upper = upper, tol = tol, ErrorL2deriv = ErrorL2deriv, 
                        Regressor = Regressor, risk = risk, neighbor = neighbor, 
                        z.comp = z.comp, stand = A, cent = z)$root, silent = TRUE)

            if(!is.numeric(b)){
                if(warn) cat("Could not determine optimal clipping bound!\n", 
                             "'radius >= maximum radius' for the given risk?\n",
                             "=> the minimum asymptotic bias (lower case) solution is returned\n")
                res <- getInfRobRegTypeIC(ErrorL2deriv = ErrorL2deriv, Regressor = Regressor, 
                                risk = asBias(), neighbor = neighbor, 
                                ErrorL2derivDistrSymm = ErrorL2derivDistrSymm, 
                                trafo = trafo, maxiter = maxiter, tol = tol, warn = warn)
                Risk <- getAsRiskRegTS(risk = risk, ErrorL2deriv = ErrorL2deriv, 
                                       Regressor = Regressor, neighbor = neighbor,
                                       clip = res$b, cent = res$a, stand = res$A, 
                                       trafo = trafo)
                res$risk <- c(Risk, res$risk)
                return(res)
            }
            z <- getInfCentRegTS(ErrorL2deriv = ErrorL2deriv, Regressor = Regressor, 
                        neighbor = neighbor, clip = b, cent = z, stand = A, z.comp = z.comp)

            if(is(Regressor, "UnivariateDistribution"))
                A <- A.old
            else
                A <- getInfStandRegTS(ErrorL2deriv = ErrorL2deriv, Regressor = Regressor,
                            neighbor = neighbor, z.comp = z.comp, clip = b, cent = z, 
                            stand = A, trafo = trafo)

            prec <- max(abs(b-b.old), abs(A-A.old), abs(z-z.old))
#            cat("current precision in IC algo:\t", prec, "\n")
            if(is(Regressor, "UnivariateDistribution") & (!z.comp)) break
            if(prec < tol) break
            if(iter > maxiter){
                cat("maximum iterations reached!\n", "achieved precision:\t", prec, "\n")
                break
            }
        }
        if(is(Regressor, "UnivariateDistribution")){
            A <- getInfStandRegTS(ErrorL2deriv = ErrorL2deriv, Regressor = Regressor,
                        neighbor = neighbor, z.comp = z.comp, clip = b, cent = z, 
                        stand = A, trafo = trafo)
            b <- as.vector(A) * b
        }
        
        a <- as.vector(A %*% z)
        info <- paste("optimally robust IC for", sQuote(class(risk)[1]))
        Risk <- getAsRiskRegTS(risk = risk, ErrorL2deriv = ErrorL2deriv, 
                               Regressor = Regressor, neighbor = neighbor,
                               clip = b, cent = a, stand = A, trafo = trafo)
        Risk <- c(Risk, list(asBias = b))

        return(list(A = A, a = a, b = b, d = NULL, risk = Risk, info = info))    
    })
setMethod("getInfRobRegTypeIC", signature(ErrorL2deriv = "RealRandVariable",
                                          Regressor = "Distribution", 
                                          risk = "asGRisk", 
                                          neighbor = "ContNeighborhood"),
    function(ErrorL2deriv, Regressor, risk, neighbor, ErrorSymm, RegSymm, 
             ErrorDistr, ErrorL2derivSymm, ErrorL2derivDistrSymm, 
             Finfo, trafo, upper, z.start, A.start, maxiter, tol, warn){
        if(is.null(z.start)) z.start <- numeric(ncol(trafo))
        if(is.null(A.start)) A.start <- trafo

        radius <- neighbor@radius
        if(identical(all.equal(radius, 0), TRUE)){
            if(warn) cat("'radius == 0' => (classical) optimal IC\n", 
                         "in sense of Cramer-Rao bound is returned\n")
            res <- getInfRobRegTypeIC(ErrorL2deriv = ErrorL2deriv, Regressor = Regressor,
                        risk = asCov(), neighbor = neighbor, ErrorDistr = ErrorDistr, 
                        Finfo = Finfo, trafo = trafo)
            Risk <- getAsRiskRegTS(risk = risk, ErrorL2deriv = ErrorL2deriv, 
                                   Regressor = Regressor, neighbor = neighbor,
                                   clip = res$b, cent = res$a, stand = res$A, 
                                   trafo = trafo)
            res$risk <- c(Risk, res$risk)
            return(res)
        }
        z <- z.start
        A <- A.start
        b <- 0
        nrvalues <- length(ErrorL2deriv)
        z.comp <- rep(TRUE, nrvalues)
        A.comp <- matrix(TRUE, ncol = nrvalues, nrow = nrvalues)
        for(i in 1:nrvalues){
            if(is(ErrorL2derivDistrSymm[[i]], "SphericalSymmetry"))
                if(ErrorL2derivDistrSymm[[i]]@SymmCenter == 0)
                    z.comp[i] <- FALSE
        }

        for(i in 1:(nrvalues-1))
            for(j in (i+1):nrvalues){
                if(is(ErrorSymm, "SphericalSymmetry")){
                    if((is(ErrorL2derivSymm[[i]], "OddSymmetric") & is(ErrorL2derivSymm[[j]], "EvenSymmetric"))
                       | (is(ErrorL2derivSymm[[j]], "OddSymmetric") & is(ErrorL2derivSymm[[i]], "EvenSymmetric")))
                        if((ErrorL2derivSymm[[i]]@SymmCenter == ErrorL2derivSymm[[j]]@SymmCenter)
                           & (ErrorL2derivSymm[[i]]@SymmCenter == ErrorSymm@SymmCenter))
                            A.comp[i,j] <- FALSE
                }
            }
        A.comp[col(A.comp) < row(A.comp)] <- A.comp[col(A.comp) > row(A.comp)] 

        iter <- 0
        repeat{
            iter <- iter + 1
            z.old <- z
            b.old <- b
            A.old <- A
            tb <- system.time(b <- try(uniroot(getInfClipRegTS, lower = .Machine$double.eps^0.75, 
                        upper = upper, tol = tol, ErrorL2deriv = ErrorL2deriv, 
                        Regressor = Regressor, risk = risk, neighbor = neighbor, 
                        ErrorDistr = ErrorDistr, stand = A, cent = z, trafo = trafo)$root, silent = FALSE))
            cat("time for b-step:\t", tb, "\n")
            if(!is.numeric(b)){
                if(warn) cat("Could not determine optimal clipping bound!\n", 
                             "'radius >= maximum radius' for the given risk?\n",
                             "=> the minimum asymptotic bias (lower case) solution is returned\n",
                             "If 'no' => Try again with modified starting values ",
                             "'z.start' and 'A.start'\n") 
                res <- getInfRobRegTypeIC(ErrorL2deriv = ErrorL2deriv, Regressor = Regressor,
                                risk = asBias(), neighbor = neighbor, ErrorSymm = ErrorSymm,
                                RegSymm = RegSymm, ErrorDistr = ErrorDistr, 
                                ErrorL2derivSymm = ErrorL2derivSymm, 
                                ErrorL2derivDistrSymm = ErrorL2derivDistrSymm, 
                                Finfo = Finfo, trafo = trafo, z.start = z.start, A.start = A.start,
                                upper = upper, maxiter = maxiter, tol = tol, warn = warn)
                Risk <- getAsRiskRegTS(risk = risk, ErrorL2deriv = ErrorL2deriv, 
                                       Regressor = Regressor, neighbor = neighbor,
                                       clip = res$b, cent = res$a, stand = res$A, 
                                       trafo = trafo)
                res$risk <- c(Risk, res$risk)
                return(res)
            }
            tz <- system.time(z <- getInfCentRegTS(ErrorL2deriv = ErrorL2deriv, Regressor = Regressor, 
                            neighbor = neighbor, ErrorDistr = ErrorDistr, stand = A, 
                            cent = z, clip = b, z.comp = z.comp))
            cat("time for z-step:\t", tz, "\n")
            tA <- system.time(A <- getInfStandRegTS(ErrorL2deriv = ErrorL2deriv, Regressor = Regressor, 
                        neighbor = neighbor, ErrorDistr = ErrorDistr, A.comp = A.comp,
                        stand = A, clip = b, cent = z, trafo = trafo))
            cat("time for A-step:\t", tA, "\n")

            prec <- max(abs(b-b.old), max(abs(A-A.old)), max(abs(z-z.old)))
            cat("current precision in IC algo:\t", prec, "\n")
            if(prec < tol) break
            if(iter > maxiter){
                cat("maximum iterations reached!\n", "achieved precision:\t", prec, "\n")
                break
            }
        }
        a <- as.vector(A %*% z)
        info <- paste("optimally robust IC for", sQuote(class(risk)[1]))
        Risk <- getAsRiskRegTS(risk = risk, ErrorL2deriv = ErrorL2deriv, 
                               Regressor = Regressor, neighbor = neighbor,
                               clip = b, cent = a, stand = A, trafo = trafo)
        Risk <- c(Risk, list(asBias = b))

        return(list(A = A, a = a, b = b, d = NULL, risk = Risk, info = info))    
    })
