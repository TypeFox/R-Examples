###############################################################################
## IC algorithm for asymptotic Anscombe risk
###############################################################################
setMethod("getInfRobIC", signature(L2deriv = "UnivariateDistribution", 
                                   risk = "asAnscombe", 
                                   neighbor = "UncondNeighborhood"),
    function(L2deriv, risk, neighbor, symm, Finfo, trafo, 
             upper = NULL, lower = NULL, maxiter, tol, warn, noLow = FALSE,
             verbose = NULL, checkBounds = TRUE, ...){

        if(missing(verbose)|| is.null(verbose))
           verbose <- getRobAStBaseOption("all.verbose")

        eff <- eff(risk)
        
        i <- 0

        maxi <- min(5,maxiter%/%4)
        toli <- min(tol*100,1e-3)

        FI0 <- trafo%*%matrix(1/Finfo)%*%t(trafo)
        normtype <- normtype(risk)
        std <- if(is(normtype(risk),"QFNorm")) 
                      QuadForm(normtype(risk)) else diag(nrow(trafo))
        FI <- sum(diag(std%*%FI0))
                
        
        
        lowBerg <- minmaxBias(L2deriv = L2deriv, neighbor = neighbor,
                   biastype = biastype(risk), symm = symm, 
                   trafo = trafo, maxiter = maxi, 
                   tol = toli, warn = warn, Finfo = Finfo)
        V <- lowBerg$risk$asCov
        trV <- sum(diag(std%*%V))
        
        
        if(FI/trV >eff){
           lowBerg$eff <- FI/trV
           return(lowBerg)
        }
        
        it.erg <- 0
        erg <- 0
        if(is.null(lower) || lower < lowBerg$risk$asBias$value) 
           { lower <- lowBerg$risk$asBias$value
             f.low <- FI/trV-eff
           } else f.low <- NULL        
        
        if(is.null(upper))
           upper <- max(4*lower,q(L2deriv)(eff^.5)*3)
  
        e.up <- 0
        while(e.up < eff){
             risk.b <- asHampel(bound = upper, biastype = biastype(risk), 
                             normtype = normtype(risk))
             upBerg <- getInfRobIC(L2deriv, risk.b, neighbor, symm, Finfo, trafo, 
                                   upper = 3*upper, lower = lower, maxiter = maxi, 
                                   tol = toli, warn, noLow = noLow,
                                   verbose = FALSE, checkBounds = FALSE) 
             trV <- upBerg$risk$trAsCov$value
             e.up <- FI/trV
             upper <- upper * 3
           } 
        upper <- upper / 3
  

        funb <- function(b0){
          risk.b <- asHampel(bound = b0, biastype = biastype(risk), 
                           normtype = normtype(risk))
          it.erg <<- it.erg + 1
          maxi <- min(5,maxiter%/%4^(1/it.erg))
          toli <- min(tol*100^(1/it.erg),1e-3)
          checkBounds <- checkBounds & it.erg>10
          erg <<- getInfRobIC(L2deriv, risk.b, neighbor, symm, Finfo, trafo, 
             upper = upper, lower = lower, maxiter = maxi, tol = toli, warn, noLow = noLow,
             verbose = verbose, checkBounds = checkBounds)
          trV <- erg$risk$trAsCov$value
          if(verbose) cat("Outer iteration:", it.erg,"  b_0=", round(b0,3), 
                          " eff=", round(FI/trV,3), "\n")  
          return(FI/trV-eff)
          }

        if(is.null(f.low)) f.low  <- fun(lower)  

        if(verbose) print(c(lower,upper, f.lower=f.low, f.upper=e.up-eff))
        
        b <- uniroot(funb, interval=c(lower,upper), f.lower=f.low, 
                     f.upper=e.up-eff,tol=tol,maxiter=maxiter)
        erg$info <- c(erg$info,
                  paste("optimally bias-robust IC for ARE", eff, " in the ideal model"))

        erg$risk$eff <- b$f.root+eff
        return(erg)
        })

###################################################################################
# multivariate solution Anscombe   --- new 24-08-10
###################################################################################

setMethod("getInfRobIC", signature(L2deriv = "RealRandVariable",
                                   risk = "asAnscombe",
                                   neighbor = "UncondNeighborhood"),
    function(L2deriv, risk, neighbor, Distr, DistrSymm, L2derivSymm,
             L2derivDistrSymm, Finfo, trafo, onesetLM = FALSE,
             z.start, A.start, upper = NULL, lower = NULL,
             OptOrIter = "iterate", maxiter, tol, warn,
             verbose = NULL, checkBounds = TRUE, ...){

        if(missing(verbose)|| is.null(verbose))
           verbose <- getRobAStBaseOption("all.verbose")

        mc <- match.call()

        eff <- eff(risk)

        ## some abbreviations / checks
        biastype <- biastype(risk)
        normtype <- normtype(risk)

        p <- nrow(trafo)
        k <- ncol(trafo)

        maxi <- min(5,maxiter%/%4)
        toli <- min(tol*100,1e-3)

        std <- if(is(normtype,"QFNorm")) QuadForm(normtype) else diag(p)

        if(! is(neighbor,"ContNeighborhood") && p>1)
           stop("Not yet implemented")

        ## non-standard norms
        FI1 <- trafo%*%solve(Finfo)
        FI0 <- FI1%*%t(trafo)
        FI <- solve(FI0)
        if(is(normtype,"InfoNorm") || is(normtype,"SelfNorm") ){
           QuadForm(normtype) <- PosSemDefSymmMatrix(FI)
           normtype(risk) <- normtype
        }
        std <- if(is(normtype,"QFNorm"))
                  QuadForm(normtype) else diag(p)
        trV.ML <- sum(diag(std%*%FI0))

        if(is.null(upper))
           upper <- sqrt(eff*max(diag(std%*%FI0)))*3
        
        
        lowBerg <- .getLowerSol(L2deriv = L2deriv, risk = risk,
                                   neighbor = neighbor, Distr = Distr,
                                   DistrSymm = DistrSymm,
                                   L2derivSymm = L2derivSymm,
                                   L2derivDistrSymm = L2derivDistrSymm,
                                   z.start = z.start, A.start = A.start,
                                   trafo = trafo, maxiter = maxiter, 
                                   tol = tol,
                                   warn = FALSE, Finfo = Finfo, 
                                   QuadForm = std, verbose = verbose)

        if(is.null(lower)||(lower< lowBerg$b))
           {lower <- lowBerg$b
#            print(lowBerg$risk$asAnscombe)
            f.low <- lowBerg$risk$asAnscombe - eff 
        } else {
             risk.b <- asHampel(bound = lower, biastype = biastype, 
                             normtype = normtype)            
             lowBerg <- getInfRobIC(L2deriv, risk.b, neighbor, 
                 Distr, DistrSymm, L2derivSymm,
                 L2derivDistrSymm, Finfo, trafo, onesetLM = onesetLM,
                 z.start, A.start, upper = upper, lower = lower,
                 OptOrIter = OptOrIter, maxiter=maxi, 
                 tol=toli, warn,
                 verbose = FALSE, checkBounds = FALSE, ...)
             trV <- lowBerg$risk$trAsCov$value
             f.low <- trV.ML/trV -eff 
        }
        
        if(f.low > 0){
           lowBerg$call <- mc 
           lowBerg$eff <- f.low + eff
           return(lowBerg)
        }

        e.up <- 0
        if(lower>=upper) upper <- lower*3
        while(e.up < eff){
           risk.b <- asHampel(bound = upper, biastype = biastype, 
                             normtype = normtype)
           upBerg <- getInfRobIC(L2deriv, risk.b, neighbor, 
                 Distr, DistrSymm, L2derivSymm,
                 L2derivDistrSymm, Finfo, trafo, onesetLM = onesetLM,
                 z.start, A.start, upper = upper, lower = lower,
                 OptOrIter = OptOrIter, maxiter=maxi, 
                 tol=toli, warn,
                 verbose = FALSE, checkBounds = FALSE, ...)
           trV <- upBerg$risk$trAsCov$value
           e.up <- trV.ML/trV
           upper <- upper * 3
        } 
        upper <- upper / 3


        erg <- 0
        it.erg <- 0
        funb <- function(b0){
          risk.b <- asHampel(bound = b0, biastype = biastype(risk), 
                             normtype = normtype(risk))
          it.erg <<- it.erg + 1
          maxi <- min(5,maxiter%/%4^(1/it.erg))
          toli <- min(tol*100^(1/it.erg),1e-3)
          chkbd <- if(it.erg<25) FALSE else checkBounds
          verbL <- if(it.erg<25) FALSE else verbose
          
          erg <<- getInfRobIC(L2deriv, risk.b, neighbor, 
             Distr, DistrSymm, L2derivSymm,
             L2derivDistrSymm, Finfo, trafo, onesetLM = onesetLM,
             z.start, A.start, upper = upper, lower = lower,
             OptOrIter = OptOrIter, maxiter = maxi, tol = toli , warn,
             verbose = verbL, checkBounds = chkbd, ...)
          trV <- erg$risk$trAsCov$value
          if(verbose) cat("Outer iteration:", it.erg,"  b_0=", round(b0,3), 
                          " eff=", round(trV.ML/trV,3), "\n")  
          return(trV.ML/trV-eff)
          }
        print(c(lower,upper, f.lower=f.low, f.upper=e.up-eff))
        b <- uniroot(funb, interval=c(lower,upper), f.lower=f.low, 
                     f.upper=e.up-eff,tol=tol,maxiter=maxiter)
        erg$info <- c(erg$info,
                  paste("optimally bias-robust IC for ARE", eff, " in the ideal model"))

        erg$risk$eff <- b$f.root+eff
        erg$call <- mc 
        return(erg)
        }
        )

