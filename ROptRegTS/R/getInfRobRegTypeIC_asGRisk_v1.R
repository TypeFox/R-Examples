###############################################################################
## get optimally robust IC for convex asymptotic risks
###############################################################################
setMethod("getInfRobRegTypeIC", signature(ErrorL2deriv = "UnivariateDistribution",
                                          Regressor = "Distribution", 
                                          risk = "asGRisk", 
                                          neighbor = "Av1CondTotalVarNeighborhood"),
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
        z <- function(x){numeric(1)}

        A <- trafo
        b <- 0
        if(is(ErrorL2derivDistrSymm, "SphericalSymmetry")) 
            z.comp <- !(ErrorL2derivDistrSymm@SymmCenter == 0)
        else
            z.comp <- TRUE
        
        if(is(Regressor, "UnivariateDistribution")){
            if(is(Regressor, "AbscontDistribution")){
                xlower <- ifelse(is.finite(q(Regressor)(0)), q(Regressor)(0), q(Regressor)(distr::TruncQuantile))
                xupper <- ifelse(is.finite(q(Regressor)(1)), q(Regressor)(1), q(Regressor)(1 - distr::TruncQuantile))
                x.vec <- seq(from = xlower, to = xupper, length = 1000)
            }else{
                if(is(Regressor, "DiscreteDistribution"))
                    x.vec <- support(Regressor) 
                else
                    x.vec <- unique(r(Regressor)(distr::RtoDPQ.e))
            }
            z.vec <- numeric(length(x.vec))
        }else{
            if(is(Regressor, "DiscreteMVDistribution"))
                x.vec <- support(Regressor)
            else
                stop("not yet implemented")
    
            z.vec <- numeric(nrow(x.vec))
        }

        iter <- 0
        repeat{
            iter <- iter + 1
            b.old <- b
            z.vec.old <- z.vec
            A.old <- A
            b <- try(uniroot(getInfClipRegTS, lower = .Machine$double.eps^0.75, 
                        upper = upper, tol = tol, ErrorL2deriv = ErrorL2deriv, 
                        Regressor = Regressor, risk = risk, neighbor = neighbor, 
                        z.comp = z.comp, stand = A, cent = z)$root, silent = FALSE)
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
            zz <- getInfCentRegTS(ErrorL2deriv = ErrorL2deriv, Regressor = Regressor, 
                        neighbor = neighbor, clip = b, cent = z, stand = A, 
                        z.comp = z.comp, x.vec = x.vec, tol.z = tol)
            z.vec <- zz$z.vec
            z <- zz$z.fct
            if(is(Regressor, "UnivariateDistribution"))
                A <- A.old
            else
                A <- getInfStandRegTS(ErrorL2deriv = ErrorL2deriv, Regressor = Regressor,
                            neighbor = neighbor, z.comp = z.comp, clip = b, cent = z, 
                            stand = A, trafo = trafo)
            prec <- max(abs(b-b.old), abs(A-A.old), abs(z.vec-z.vec.old))
            if(z.comp)
                cat("current precision in IC algo:\t", prec, "\n")
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
            a.fct <- function(x){ A * z(x) }
            body(a.fct) <- substitute({ A * z(x) }, list(z = z))
            a <- RealRandVariable(Map = list(a.fct), Domain = img(Regressor))
        }else{
            a <- RealRandVariable(Map = list(z), Domain = img(Regressor))
        }

        info <- paste("optimally robust IC for", sQuote(class(risk)[1]))
        Risk <- getAsRiskRegTS(risk = risk, ErrorL2deriv = ErrorL2deriv, 
                               Regressor = Regressor, neighbor = neighbor,
                               clip = b, cent = a, stand = A, trafo = trafo)
        Risk <- c(Risk, list(asBias = b))

        return(list(A = A, a = a, b = b, d = NULL, risk = Risk, info = info))    
    })
