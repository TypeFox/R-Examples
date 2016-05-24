###############################################################################
## get optimally robust IC for asymptotic under-/overshoot risk
###############################################################################
setMethod("getInfRobIC", signature(L2deriv = "UnivariateDistribution", 
                                   risk = "asUnOvShoot", 
                                   neighbor = "UncondNeighborhood"),
    function(L2deriv, risk, neighbor, symm, Finfo, trafo, 
            upper, maxiter, tol, warn){
        radius <- neighbor@radius
        if(identical(all.equal(radius, 0), TRUE)){
            if(warn) cat("'radius == 0' => (classical) optimal IC\n", 
                         "in sense of Cramer-Rao bound is returned\n")
            res <- getInfRobIC(L2deriv = L2deriv, risk = asCov(), 
                        neighbor = TotalVarNeighborhood(radius = neighbor@radius), 
                        Finfo = Finfo, trafo = trafo)
            Risk <- getAsRisk(risk = risk, L2deriv = L2deriv, neighbor = neighbor, 
                              clip = res$b, cent = res$a, stand = res$A, trafo = trafo)
            res$risk <- c(Risk, res$risk)
            return(res)
        }
        
        bound <- risk@width*(-m1df(L2deriv, 0))
        if(is(neighbor, "ContNeighborhood")){
            if(identical(all.equal(radius, 2*bound), TRUE)){
                zi <- sign(as.vector(trafo))
                A <- as.matrix(zi)
                b <- zi*as.vector(trafo)*2*risk@width/radius
                p0 <- p(L2deriv)(0)
                if(is(L2deriv, "AbscontDistribution"))
                    ws0 <- 0
                else
                    ws0 <- d(L2deriv)(0)

                if(zi == 1)
                    a <- -b*(1-p0)/(1-ws0)
                else
                    a <- b*(p0-ws0)/(1-ws0)
            
                info <- paste("optimally robust IC for", sQuote(class(risk)[1]))
                Risk <- list(asUnOvShoot = 0.5)

                return(list(A = A, a = a, b = b, d = 1, risk = Risk, info = info))                
            }
            if(radius > 2*bound)
                stop("boundedness condition is violated!")
        }

        if(is(neighbor, "TotalVarNeighborhood")){
            if(identical(all.equal(radius, bound), TRUE)){
                zi <- sign(as.vector(trafo))
                A <- as.matrix(zi)
                b <- zi*as.vector(trafo)*risk@width/radius
                p0 <- p(L2deriv)(0)
                if(is(L2deriv, "AbscontDistribution"))
                    ws0 <- 0
                else
                    ws0 <- d(L2deriv)(0)

                if(zi == 1)
                    a <- -b*(1-p0)/(1-ws0)
                else
                    a <- b*(p0-ws0)/(1-ws0)
            
                info <- paste("optimally robust IC for", sQuote(class(risk)[1]))
                Risk <- list(asUnOvShoot = 0.5)

                return(list(A = A, a = a, b = b, d = 1, risk = Risk, info = info))                
            }
            if(radius > bound)
                stop("boundedness condition is violated!")
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
            z <- getInfCent(L2deriv = L2deriv, neighbor = TotalVarNeighborhood(radius = neighbor@radius), 
                        clip = c0, cent = z, symm = S, trafo = trafo, tol.z = tol)
#            cat("c0:\t", c0, "c0.old:\t", c0.old, "z:\t", z, "z.old:\t", z.old, "\n")
            if(S) break
            if(max(abs(z - z.old), abs(c0-c0.old)) < tol) break
            if(iter > maxiter){
                cat("maximum iterations reached!\n", "achieved precision:\t", abs(c0 - c0.old), "\n")
                break
            }
        }
        info <- paste("optimally robust IC for", sQuote(class(risk)[1]))
        A <- getInfStand(L2deriv = L2deriv, neighbor = TotalVarNeighborhood(radius = neighbor@radius), 
                    clip = c0, cent = z, trafo = trafo)
        a <- as.vector(A)*z
        b <- abs(as.vector(A))*c0
        Risk <- getAsRisk(risk = risk, L2deriv = L2deriv, neighbor = neighbor, 
                          clip = b, cent = a, stand = A, trafo = trafo)

        return(list(A = A, a = a, b = b, d = NULL, risk = Risk, info = info))    
    })
