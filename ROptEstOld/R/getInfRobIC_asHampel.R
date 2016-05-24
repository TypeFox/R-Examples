###############################################################################
## IC algorithm for asymptotic Hampel risk
###############################################################################
setMethod("getInfRobIC", signature(L2deriv = "UnivariateDistribution", 
                                   risk = "asHampel", 
                                   neighbor = "UncondNeighborhood"),
    function(L2deriv, risk, neighbor, symm, Finfo, trafo, 
             upper, maxiter, tol, warn){
        A <- trafo / E(L2deriv, function(x){x^2})
        b <- risk@bound
        bmax <- abs(as.vector(A))*max(abs(q(L2deriv)(0)), q(L2deriv)(1))
        if(b >= bmax){
            if(warn) cat("'b >= maximum asymptotic bias' => (classical) optimal IC\n", 
                         "in sense of Cramer-Rao bound is returned\n")
            res <- getInfRobIC(L2deriv = L2deriv, risk = asCov(), 
                        neighbor = neighbor, Finfo = Finfo, trafo = trafo)
            return(res)
        }

        bmin <- getAsRisk(risk = asBias(), L2deriv = L2deriv, neighbor = neighbor, 
                          trafo = trafo)$asBias
        if(b <= bmin){
            if(warn) cat("'b <= minimum asymptotic bias'\n",
                         "=> the minimum asymptotic bias (lower case) solution is returned\n")
            res <- getInfRobIC(L2deriv = L2deriv, risk = asBias(), 
                            neighbor = neighbor, symm = symm, 
                            trafo = trafo, maxiter = maxiter, tol = tol)
            Risk <- list(asMSE = res$risk$asCov + neighbor@radius^2*bmin^2)
            res$risk <- c(Risk, res$risk)
            return(res)
        }
        c0 <- b/as.vector(A)
        if(is(symm, "SphericalSymmetry")) 
            S <- symm@SymmCenter == 0
        else
            S <- FALSE
        z <- getInfCent(L2deriv = L2deriv, neighbor = neighbor, 
                    clip = c0, cent = 0, trafo = trafo, tol.z = tol, symm = S)
        iter <- 0
        repeat{
            iter <- iter + 1
            A.old <- A
            z.old <- z
            A <- getInfStand(L2deriv = L2deriv, neighbor = neighbor, 
                        clip = c0, cent = z, trafo = trafo)
            c0 <- b/as.vector(A)
            z <- getInfCent(L2deriv = L2deriv, neighbor = neighbor, 
                        clip = c0, cent = z, trafo = trafo, tol.z = tol, symm = S)
            if(max(abs(as.vector(A-A.old)), abs(z-z.old)) < tol) break
            if(iter > maxiter){
                cat("maximum iterations reached!\n", "achieved precision:\t", 
                    max(abs(as.vector(A-A.old)), abs(z-z.old)), "\n")
                break
            }
        }
        info <- paste("optimally robust IC for 'asHampel' with bound =", round(b,3))
        a <- as.vector(A)*z
        Cov <- getAsRisk(risk = asCov(), L2deriv = L2deriv, neighbor = neighbor, 
                         clip = b, cent = a, stand = A)$asCov
        Risk <- list(asCov = Cov, asBias = b, asMSE = Cov + neighbor@radius^2*b^2)

        return(list(A = A, a = a, b = b, d = NULL, risk = Risk, info = info))
    })
setMethod("getInfRobIC", signature(L2deriv = "RealRandVariable", 
                                   risk = "asHampel", 
                                   neighbor = "ContNeighborhood"),
    function(L2deriv, risk, neighbor, Distr, DistrSymm, L2derivSymm,
             L2derivDistrSymm, Finfo, trafo, z.start, A.start, upper, maxiter, tol, warn){
        if(is.null(z.start)) z.start <- numeric(ncol(trafo))
        if(is.null(A.start)) A.start <- trafo

        ClassIC <- trafo %*% solve(Finfo) %*% L2deriv
        lower <- q(Distr)(getdistrOption("TruncQuantile"))
        upper <- q(Distr)(1-getdistrOption("TruncQuantile"))
        x <- seq(from = lower, to = upper, by = 0.01)
        bmax <- evalRandVar(ClassIC, as.matrix(x))^2
        bmax <- sqrt(max(colSums(bmax)))
        b <- risk@bound
        cat("numerical approximation of maximal bound:\t", bmax, "\n")
        if(b >= bmax){
            if(warn) cat("'b >= maximum asymptotic bias' => (classical) optimal IC\n", 
                         "in sense of Cramer-Rao bound is returned\n")
            res <- getInfRobIC(L2deriv = L2deriv, risk = asCov(), neighbor = neighbor, 
                               Distr = Distr, Finfo = Finfo, trafo = trafo)
            return(res)
        }
        bmin <- getAsRisk(risk = asBias(), L2deriv = L2deriv, neighbor = neighbor, 
                          Distr = Distr, L2derivDistrSymm = L2derivDistrSymm, 
                          trafo = trafo, z.start = z.start, A.start = A.start, 
                          maxiter = maxiter, tol = tol)$asBias
        cat("minimal bound:\t", bmin, "\n")
        if(b <= bmin){
            if(warn) cat("'b <= minimum asymptotic bias'\n",
                         "=> the minimum asymptotic bias (lower case) solution is returned\n")
            res <- getInfRobIC(L2deriv = L2deriv, risk = asBias(), neighbor = neighbor, 
                               Distr = Distr, DistrSymm = DistrSymm, L2derivSymm = L2derivSymm,
                               L2derivDistrSymm = L2derivDistrSymm, z.start = z.start, 
                               A.start = A.start, trafo = trafo, maxiter = maxiter, 
                               tol = tol, warn = warn)
            return(res)
        }
        nrvalues <- length(L2deriv)
        z.comp <- rep(TRUE, nrvalues)
        A.comp <- matrix(1, ncol = nrvalues, nrow = nrvalues)
        for(i in 1:nrvalues){
            if(is(L2derivDistrSymm[[i]], "SphericalSymmetry"))
                if(L2derivDistrSymm[[i]]@SymmCenter == 0){
                    z.comp[i] <- FALSE
                    A.comp[i,i] <- 2
                }
        }
        for(i in 1:(nrvalues-1))
            for(j in (i+1):nrvalues){
                if(z.comp[i] | z.comp[j]){ 
                    if(is(DistrSymm, "SphericalSymmetry")){
                        if((is(L2derivSymm[[i]], "OddSymmetric") & is(L2derivSymm[[j]], "EvenSymmetric"))
                           | (is(L2derivSymm[[j]], "OddSymmetric") & is(L2derivSymm[[i]], "EvenSymmetric")))
                            if((L2derivSymm[[i]]@SymmCenter == L2derivSymm[[j]]@SymmCenter)
                               & (L2derivSymm[[i]]@SymmCenter == DistrSymm@SymmCenter))
                                A.comp[i,j] <- 0
                    }else{
                        A.comp[i,j] <- 2
                    }
                }
            }
        A.comp[col(A.comp) < row(A.comp)] <- A.comp[col(A.comp) > row(A.comp)]

        z <- z.start
        A <- A.start
        iter <- 0
        repeat{
            iter <- iter + 1
            z.old <- z
            A.old <- A
            z <- getInfCent(L2deriv = L2deriv, neighbor = neighbor, Distr = Distr, 
                            z.comp = z.comp, stand = A, cent = z, clip = b)
            A <- getInfStand(L2deriv = L2deriv, neighbor = neighbor, clip = b, cent = z, 
                             A.comp = A.comp, trafo = trafo, Distr = Distr, stand = A)
            prec <- max(max(abs(A-A.old)), max(abs(z-z.old)))
            cat("current precision in IC algo:\t", prec, "\n")
            if(prec < tol) break
            if(iter > maxiter){
                cat("maximum iterations reached!\n", "achieved precision:\t", prec, "\n")
                break
            }
        }
        info <- paste("optimally robust IC for 'asHampel' with bound =", round(b,3))
        a <- as.vector(A %*% z)
        Cov <- getAsRisk(risk = asCov(), L2deriv = L2deriv, neighbor = neighbor, 
                           Distr = Distr, clip = b, cent = a, stand = A)$asCov
        Risk <- list(asCov = Cov, asBias = b, asMSE = sum(diag(Cov)) + neighbor@radius^2*b^2)

        return(list(A = A, a = a, b = b, d = NULL, risk = Risk, info = info))
    })
