###############################################################################
## get optimally robust IC for convex asymptotic risks
###############################################################################
setMethod("getInfRobIC", signature(L2deriv = "UnivariateDistribution", 
                                   risk = "asGRisk", 
                                   neighbor = "UncondNeighborhood"),
    function(L2deriv, risk, neighbor, symm, Finfo, trafo, 
            upper, maxiter, tol, warn){
        radius <- neighbor@radius
        if(identical(all.equal(radius, 0), TRUE)){
            if(warn) cat("'radius == 0' => (classical) optimal IC\n", 
                         "in sense of Cramer-Rao bound is returned\n")
            res <- getInfRobIC(L2deriv = L2deriv, risk = asCov(), 
                        neighbor = neighbor, Finfo = Finfo, trafo = trafo)
            Risk <- getAsRisk(risk = risk, L2deriv = L2deriv, neighbor = neighbor, 
                              clip = res$b, cent = res$a, stand = res$A, trafo = trafo)
            res$risk <- c(Risk, res$risk)
            return(res)
        }
        z <- 0
        c0 <- 0
        iter <- 0
        if(is(symm, "SphericalSymmetry")) 
            S <- symm@SymmCenter == 0
        else
            S <- FALSE
                
        repeat{
            iter <- iter + 1
            z.old <- z
            c0.old <- c0
            c0 <- try(uniroot(getInfClip, lower = .Machine$double.eps^0.75, 
                        upper = upper, tol = tol, L2deriv = L2deriv, risk = risk, 
                        neighbor = neighbor, cent = z, symm = S, 
                        trafo = trafo)$root, silent = TRUE)
            if(!is.numeric(c0)){
                if(warn) cat("The IC algorithm did not converge!\n", 
                             "'radius >= maximum radius' for the given risk?\n",
                             "=> the minimum asymptotic bias (lower case) solution is returned\n")
                res <- getInfRobIC(L2deriv = L2deriv, risk = asBias(), 
                                neighbor = neighbor, Finfo = Finfo, 
                                symm = symm, trafo = trafo, upper = upper, 
                                maxiter = maxiter, tol = tol, warn = warn)
                Risk <- getAsRisk(risk = risk, L2deriv = L2deriv, neighbor = neighbor, 
                                  clip = res$b, cent = res$a, stand = res$A, trafo = trafo)
                res$risk <- c(Risk, res$risk)
                return(res)
            }
            z <- getInfCent(L2deriv = L2deriv, neighbor = neighbor, clip = c0, 
                        cent = z, symm = S, trafo = trafo, tol.z = tol)
#            cat("c0:\t", c0, "c0.old:\t", c0.old, "z:\t", z, "z.old:\t", z.old, "\n")
            if(S) break
            if(max(abs(z - z.old), abs(c0-c0.old)) < tol) break
            if(iter > maxiter){
                cat("maximum iterations reached!\n", "achieved precision:\t", abs(c0 - c0.old), "\n")
                break
            }
        }
        info <- paste("optimally robust IC for", sQuote(class(risk)[1]))
        A <- getInfStand(L2deriv = L2deriv, neighbor = neighbor, 
                    clip = c0, cent = z, trafo = trafo)
        a <- as.vector(A)*z
        b <- abs(as.vector(A))*c0
        Risk <- getAsRisk(risk = risk, L2deriv = L2deriv, neighbor = neighbor, 
                          clip = b, cent = a, stand = A, trafo = trafo)
        Risk <- c(Risk, list(asBias = b))

        return(list(A = A, a = a, b = b, d = NULL, risk = Risk, info = info))    
    })
setMethod("getInfRobIC", signature(L2deriv = "RealRandVariable", 
                                   risk = "asGRisk", 
                                   neighbor = "ContNeighborhood"),
    function(L2deriv, risk, neighbor, Distr, DistrSymm, L2derivSymm, 
             L2derivDistrSymm, Finfo, trafo, z.start, A.start, upper, 
             maxiter, tol, warn){
        if(is.null(z.start)) z.start <- numeric(ncol(trafo))
        if(is.null(A.start)) A.start <- trafo

        radius <- neighbor@radius
        if(identical(all.equal(radius, 0), TRUE)){
            if(warn) cat("'radius == 0' => (classical) optimal IC\n", 
                         "in sense of Cramer-Rao bound is returned\n")
            res <- getInfRobIC(L2deriv = L2deriv, risk = asCov(), neighbor = neighbor, 
                               Distr = Distr, Finfo = Finfo, trafo = trafo)
            Risk <- getAsRisk(risk = risk, L2deriv = L2deriv, neighbor = neighbor, 
                              clip = res$b, cent = res$a, stand = res$A, trafo = trafo)
            res$risk <- c(Risk, res$risk)
            return(res)
        }
        nrvalues <- length(L2deriv)
        z.comp <- rep(TRUE, nrvalues)
        A.comp <- matrix(TRUE, ncol = nrvalues, nrow = nrvalues)
        for(i in 1:nrvalues){
            if(is(L2derivDistrSymm[[i]], "SphericalSymmetry"))
                if(L2derivDistrSymm[[i]]@SymmCenter == 0)
                    z.comp[i] <- FALSE
        }
        for(i in 1:(nrvalues-1))
            for(j in (i+1):nrvalues){
                if(is(DistrSymm, "SphericalSymmetry")){
                    if((is(L2derivSymm[[i]], "OddSymmetric") & is(L2derivSymm[[j]], "EvenSymmetric"))
                       | (is(L2derivSymm[[j]], "OddSymmetric") & is(L2derivSymm[[i]], "EvenSymmetric")))
                        if((L2derivSymm[[i]]@SymmCenter == L2derivSymm[[j]]@SymmCenter)
                           & (L2derivSymm[[i]]@SymmCenter == DistrSymm@SymmCenter))
                            A.comp[i,j] <- FALSE
                }
            }
        A.comp[col(A.comp) < row(A.comp)] <- A.comp[col(A.comp) > row(A.comp)]

        z <- z.start
        A <- A.start
        b <- 0
        iter <- 0
        repeat{
            iter <- iter + 1
            z.old <- z
            b.old <- b
            A.old <- A
            b <- try(uniroot(getInfClip, lower = .Machine$double.eps^0.75, 
                        upper = upper, tol = tol, L2deriv = L2deriv, risk = risk, 
                        Distr = Distr, neighbor = neighbor, stand = A, cent = z, 
                        trafo = trafo)$root, silent = TRUE)
            if(!is.numeric(b)){
                if(warn) cat("Could not determine optimal clipping bound!\n", 
                             "'radius >= maximum radius' for the given risk?\n",
                             "=> the minimum asymptotic bias (lower case) solution is returned\n",
                             "If 'no' => Try again with modified starting values ",
                             "'z.start' and 'A.start'\n")
                res <- getInfRobIC(L2deriv = L2deriv, risk = asBias(), 
                                neighbor = neighbor, Distr = Distr, DistrSymm = DistrSymm,
                                L2derivSymm = L2derivSymm, L2derivDistrSymm = L2derivDistrSymm,
                                z.start = z.start, A.start = A.start, trafo = trafo, 
                                upper = upper, maxiter = maxiter, tol = tol, warn = warn)
                Risk <- getAsRisk(risk = risk, L2deriv = L2deriv, neighbor = neighbor, 
                                  clip = res$b, cent = res$a, stand = res$A, trafo = trafo)
                res$risk <- c(Risk, res$risk)
                return(res)
            }
            z <- getInfCent(L2deriv = L2deriv, neighbor = neighbor, Distr = Distr, 
                            z.comp = z.comp, stand = A, cent = z, clip = b)
            A <- getInfStand(L2deriv = L2deriv, neighbor = neighbor, 
                        Distr = Distr, A.comp = A.comp, stand = A, clip = b, cent = z, 
                        trafo = trafo)
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
        Risk <- getAsRisk(risk = risk, L2deriv = L2deriv, neighbor = neighbor, 
                          clip = b, cent = a, stand = A, trafo = trafo)
        Risk <- c(Risk, list(asBias = b))

        return(list(A = A, a = a, b = b, d = NULL, risk = Risk, info = info))    
    })
