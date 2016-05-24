###############################################################################
## get optimally robust IC for convex asymptotic risks
###############################################################################
setMethod("getInfRobRegTypeIC", signature(ErrorL2deriv = "UnivariateDistribution",
                                          Regressor = "Distribution", 
                                          risk = "asGRisk", 
                                          neighbor = "Av2CondContNeighborhood"),
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
                                   clip = res$b, cent = res$z, stand = res$A, 
                                   trafo = trafo)
            res$risk <- c(Risk, res$risk)
            return(res)
        }
        z <- 0
        A <- 1
        c0 <- 0
        if(is(ErrorL2derivDistrSymm, "SphericalSymmetry")) 
            z.comp <- !(ErrorL2derivDistrSymm@SymmCenter == 0)
        else
            z.comp <- TRUE

        iter <- 0
        repeat{
            iter <- iter + 1
            c0.old <- c0
            z.old <- z
            A.old <- A
            c0 <- try(uniroot(getInfClipRegTS, lower = .Machine$double.eps^0.75, 
                        upper = upper, tol = tol, ErrorL2deriv = ErrorL2deriv, 
                        Regressor = Regressor, risk = risk, neighbor = neighbor, 
                        z.comp = z.comp, stand = A, cent = z)$root, silent = FALSE)

            if(!is.numeric(c0)){
                if(warn) cat("Could not determine optimal clipping bound!\n", 
                             "'radius >= maximum radius' for the given risk?\n",
                             "=> the minimum asymptotic bias (lower case) solution is returned\n")
                res <- getInfRobRegTypeIC(ErrorL2deriv = ErrorL2deriv, Regressor = Regressor, 
                                risk = asBias(), neighbor = neighbor, 
                                ErrorL2derivDistrSymm = ErrorL2derivDistrSymm, 
                                trafo = trafo, maxiter = maxiter, tol = tol, warn = warn)
                Risk <- getAsRiskRegTS(risk = risk, ErrorL2deriv = ErrorL2deriv, 
                                       Regressor = Regressor, neighbor = neighbor,
                                       clip = res$b, cent = res$z, stand = res$A, 
                                       trafo = trafo)
                res$risk <- c(Risk, res$risk)
                return(res)
            }
            z <- getInfCentRegTS(ErrorL2deriv = ErrorL2deriv, Regressor = Regressor, 
                        neighbor = neighbor, clip = c0, cent = z, stand = A, 
                        z.comp = z.comp, tol.z = tol)

            prec <- max(abs(c0-c0.old), abs(z-z.old))
#            cat("current precision in IC algo:\t", prec, "\n")
            if(!z.comp) break
            if(prec < tol) break
            if(iter > maxiter){
                cat("maximum iterations reached!\n", "achieved precision:\t", prec, "\n")
                break
            }
        }
        A <- getInfStandRegTS(ErrorL2deriv = ErrorL2deriv, Regressor = Regressor,
                    neighbor = neighbor, z.comp = z.comp, clip = c0, cent = z, 
                    stand = A, trafo = trafo)
        b <- c0*A*sqrt(sum(diag(solve(E(Regressor, fun = function(x){ x %*% t(x) })))))
        
        info <- paste("optimally robust IC for", sQuote(class(risk)[1]))
        Risk <- getAsRiskRegTS(risk = risk, ErrorL2deriv = ErrorL2deriv, 
                               Regressor = Regressor, neighbor = neighbor,
                               clip = c0, cent = z, stand = A, trafo = trafo)
        Risk <- c(Risk, list(asBias = b))

        return(list(A = A, z = z, b = b, d = NULL, risk = Risk, info = info))    
    })
