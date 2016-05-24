###############################################################################
## IC algorithm for asymptotic Hampel risk
###############################################################################
setMethod("getInfRobIC", signature(L2deriv = "UnivariateDistribution", 
                                   risk = "asHampel", 
                                   neighbor = "UncondNeighborhood"),
    function(L2deriv, risk, neighbor, symm, Finfo, trafo, 
             upper = NULL, lower = NULL, maxiter, tol, warn, noLow = FALSE,
             verbose = NULL, checkBounds = TRUE, ...){

        if(missing(verbose)|| is.null(verbose))
           verbose <- getRobAStBaseOption("all.verbose")

        biastype <- biastype(risk)
        normtype <- normtype(risk)

        A <- trafo / E(L2deriv, function(x){x^2})
        b <- risk@bound

        if(checkBounds){
        bmax <- abs(as.vector(A))*max(abs(q(L2deriv)(0)), q(L2deriv)(1))
        if(b >= bmax){
            if(warn) cat("'b >= maximum asymptotic bias' => (classical) optimal IC\n", 
                         "in sense of Cramer-Rao bound is returned\n")
            res <- getInfRobIC(L2deriv = L2deriv, risk = asCov(), 
                               neighbor = neighbor, Finfo = Finfo, trafo = trafo,
                               verbose = verbose)
            res <- c(res, list(biastype = biastype, normtype = NormType()))
            Cov <- res$risk$asCov
            r <- neighbor@radius
            res$risk$asBias <- list(value = b, biastype = biastype, 
                                   normtype = normtype, 
                                   neighbortype = class(neighbor))
            res$risk$asMSE <- list(value = Cov + r^2*b^2, 
                                   r = r,
                                   at = neighbor)
            return(res)
        }

        if(!noLow){
            res <- getInfRobIC(L2deriv = L2deriv, risk = asBias(biastype = biastype), 
                               neighbor = neighbor, symm = symm,  
                               trafo = trafo, maxiter = maxiter, tol = tol, Finfo = Finfo,
                               warn = warn, verbose = verbose)
            bmin <- res$b
            cat("minimal bound:\t", bmin, "\n")
            }else{ 
                bmin <- b/2
            }

        if(b <= bmin){
            if(warn) cat("'b <= minimum asymptotic bias'\n",
                         "=> the minimum asymptotic bias (lower case) solution is returned\n")
            Risk <- list(asMSE = res$risk$asCov + neighbor@radius^2*bmin^2)
            res$risk <- c(Risk, res$risk)
            return(res)
        }
#        bmin <- getAsRisk(risk = asBias(biastype = biastype, normtype = normtype), 
#                          L2deriv = L2deriv, neighbor = neighbor, 
#                          biastype = biastype, trafo = trafo, Finfo = Finfo,
#                          warn = warn)$asBias
#        if(b <= bmin){
#            if(warn) cat("'b <= minimum asymptotic bias'\n",
#                         "=> the minimum asymptotic bias (lower case) solution is returned\n")
#            res <- getInfRobIC(L2deriv = L2deriv, risk = asBias(biastype = biastype), 
#                            neighbor = neighbor, symm = symm,  
#                            trafo = trafo, maxiter = maxiter, tol = tol, Finfo = Finfo,
#                            warn = warn)
#            Risk <- list(asMSE = res$risk$asCov + neighbor@radius^2*bmin^2)
#            res$risk <- c(Risk, res$risk)
#            return(res)
#        }
        }
        c0 <- b/as.vector(A)
        if(is(symm, "SphericalSymmetry")) 
            S <- symm@SymmCenter == 0
        else
            S <- FALSE
        z <- getInfCent(L2deriv = L2deriv, neighbor = neighbor,  
                        biastype = biastype, clip = c0, cent = 0, 
                        trafo = trafo, tol.z = tol, symm = S)
        iter <- 0
        repeat{
            iter <- iter + 1
            A.old <- A
            z.old <- z
            A <- getInfStand(L2deriv = L2deriv, neighbor = neighbor, 
                         biastype = biastype,
                         clip = c0, cent = z, trafo = trafo)
            c0 <- b/as.vector(A)
            z <- getInfCent(L2deriv = L2deriv, neighbor = neighbor, 
                         biastype = biastype,
                         clip = c0, cent = z, trafo = trafo, tol.z = tol, symm = S)
            if(max(abs(as.vector(A-A.old)), abs(z-z.old)) < tol) break
            if(verbose && iter%%5==1){
               cat("current precision in IC algo:\t",
                    max(abs(as.vector(A-A.old)), abs(z-z.old)), "\n")
                    print(round(c(A=A,z=z),3))
            }
            if(iter > maxiter){
                cat("maximum iterations reached!\n", "achieved precision:\t", 
                    max(abs(as.vector(A-A.old)), abs(z-z.old)), "\n")
                break
            }
        }
        info <- paste("optimally robust IC for 'asHampel' with bound =", round(b,3))
        a <- as.vector(A)*z
        Cov <- getInfV(L2deriv = L2deriv, neighbor = neighbor, 
                       biastype = biastype, clip = c0, cent = z, stand = A)

        # getAsRisk(risk = asHampel(), L2deriv = L2deriv, neighbor = neighbor, 
        #          biastype = biastype, clip = b, cent = a, stand = A)$asCov

        r <- neighbor@radius
        Risk <- list(asCov = Cov,
                     asBias = list(value = b, biastype = biastype, 
                                   normtype = normtype, 
                                   neighbortype = class(neighbor)), 
                     trAsCov = list(value = Cov, normtype = normtype),
                     asMSE = list(value = Cov + r^2*b^2, 
                                  r = r,
                                  at = neighbor))

        if(is(neighbor,"ContNeighborhood")){
            w <- new("HampelWeight")
            clip(w) <- b
            cent(w) <- as.numeric(z)
            stand(w) <- A
        }else if(is(neighbor,"TotalVarNeighborhood")){
            w <- new("BdStWeight")
            clip(w) <- c(0,b)+a
            stand(w) <- A
        } 
        weight(w) <- getweight(w, neighbor = neighbor, biastype = biastype, 
                               normW = NormType())
        return(list(A = A, a = a, b = b, d = NULL, risk = Risk, info = info, 
                    w = w, biastype = biastype, normtype = NormType()))
    })

###################################################################################
# multivariate solution Hampel   --- new 10-08-09
###################################################################################

setMethod("getInfRobIC", signature(L2deriv = "RealRandVariable",
                                   risk = "asHampel",
                                   neighbor = "UncondNeighborhood"),
    function(L2deriv, risk, neighbor, Distr, DistrSymm, L2derivSymm,
             L2derivDistrSymm, Finfo, trafo, onesetLM = FALSE,
             z.start, A.start, upper = NULL, lower = NULL,
             OptOrIter = "iterate", maxiter, tol, warn,
             verbose = NULL, checkBounds = TRUE, ...,
             .withEvalAsVar = TRUE){

        if(missing(verbose)|| is.null(verbose))
           verbose <- getRobAStBaseOption("all.verbose")

        mc <- match.call()

        ## some abbreviations / checks
        biastype <- biastype(risk)
        normtype <- normtype(risk)
        b <- risk@bound

        p <- nrow(trafo)
        k <- ncol(trafo)

        if(! is(neighbor,"ContNeighborhood") && p>1)
           stop("Not yet implemented")

        ## non-standard norms
        FI <- solve(trafo%*%solve(Finfo)%*%t(trafo))
        if(is(normtype,"InfoNorm") || is(normtype,"SelfNorm") ){
           QuadForm(normtype) <- PosSemDefSymmMatrix(FI)
           normtype(risk) <- normtype
        }
        std <- if(is(normtype,"QFNorm"))
                  QuadForm(normtype) else diag(p)

        ## starting values
        if(is.null(z.start)) z.start <- numeric(k)
        if(is.null(A.start)) A.start <- trafo%*%solve(Finfo)
        a.start <- as.numeric(A.start %*% z.start)

        ## initialize
        if(is(neighbor,"ContNeighborhood")){
            w <- new("HampelWeight")
        }else{ if (is(neighbor,"TotalVarNeighborhood"))
            w <- new("BdStWeight")
        }


        ## check whether given b is in (bmin, bmax)
        if(checkBounds){
            chk <- .checkUpLow(L2deriv = L2deriv,  b = b, risk = risk,
                               neighbor = neighbor, biastype = biastype,
                               normtype = normtype, Distr = Distr, Finfo = Finfo,
                               DistrSymm = DistrSymm, L2derivSymm = L2derivSymm,
                               L2derivDistrSymm = L2derivDistrSymm,
                               z.start = z.start, A.start = A.start,
                               trafo = trafo, maxiter = maxiter, tol = tol,
                               QuadForm = std, verbose = verbose,
                               nrvalpts = 5000, warn = warn)
            if(chk$up || chk$low) return(chk$res)
        }

        ## determine which entries must be computed
        # by default everything
        z.comp <- rep(TRUE,k)
        A.comp <- matrix(rep(TRUE,k*k),nrow=k)

        # otherwise if trafo == unitMatrix may use symmetry info
        if(.isUnitMatrix(trafo)){
            comp <- .getComp(L2deriv, DistrSymm, L2derivSymm, L2derivDistrSymm)
            z.comp <- comp$"z.comp"
            A.comp <- comp$"A.comp"
        }

        ## selection of the algorithm
        pM <- pmatch(tolower(OptOrIter),c("optimize","iterate"))
        OptOrIter <- pM
        if (is.na(pM)) OptOrIter <- 1

        ## determining A,a with either optimization of iteration:
        if(OptOrIter == 1)
           erg <- getLagrangeMultByOptim(b = b, L2deriv = L2deriv, risk = risk,
                      FI = Finfo, trafo = trafo, neighbor = neighbor,
                      biastype = biastype, normtype = normtype, Distr = Distr,
                      a.start = a.start, z.start = z.start, A.start = A.start,
                      w.start = w, std = std,
                      z.comp = z.comp, A.comp = A.comp, maxiter = maxiter,
                      tol = tol, verbose = verbose, ...)
        else{
           erg <- getLagrangeMultByIter(b = b, L2deriv = L2deriv, risk = risk,
                      trafo = trafo, neighbor = neighbor, biastype = biastype,
                      normtype = normtype, Distr = Distr,
                      a.start = a.start, z.start = z.start, A.start = A.start,
                      w.start = w,
                      std = std, z.comp = z.comp,
                      A.comp = A.comp, maxiter = maxiter, tol = tol,
                      verbose = verbose)
        }

        ## read out solution
        w <- erg$w
        A <- erg$A
        a <- erg$a
        z <- erg$z
        biastype <- erg$biastype
        normtype.old <- erg$normtype.old
        normtype <- erg$normtype
        risk <- erg$risk
        iter <- erg$iter
        prec <- erg$prec
        OptIterCall <- erg$call
        std <- erg$std

        ### issue some diagnostics if wanted
        if(verbose){
           cat("Iterations needed:",iter,"\n")
           cat("Precision achieved:", prec,"\n")
        }

        ## shall Lagrange-Multipliers inside weight and outside coincide
        if (onesetLM){
            if(is(neighbor,"ContNeighborhood"))
               cent(w) <- as.numeric(z)
            if(is(neighbor,"TotalVarNeighborhood"))
               clip(w) <- c(0,b)+a
            stand(w) <- A

            weight(w) <- getweight(w, neighbor = neighbor, biastype = biastype,
                                   normW = normtype)
        }
        else normtype <- normtype.old


        ### determine Covariance of pIC
        Cov <- substitute(do.call(getInfV, args = list(L2deriv = L2deriv0,
                          neighbor = neighbor0, biastype = biastype0,
                          Distr = Distr0, V.comp = A.comp0, cent = a0,
                          stand = A0, w = w0)), list(L2deriv0 = L2deriv,
                          neighbor0 = neighbor, biastype0 = biastype,
                          Distr0 = Distr, A.comp0 = A.comp, a0 = a,
                          A0 = A, w0 = w))

        rifct <- function(std0, Cov0, rad0, b0){
                     sum(diag(std0%*%eval(Cov0))) + rad0^2 * b0^2}

        asMSE.0 <- substitute(do.call(ri.fct, args=list(std0=std1, Cov0=Cov1,
                                    rad0 = rad1, b0=b1)), list(ri.fct = rifct,
                                    std1=std, Cov1=Cov, rad1=neighbor@radius,
                                    b1=b))


        ### add some further informations for the pIC-slots info and risk
        info <- paste("optimally robust IC for 'asHampel' with bound =", round(b,3))

        trAsCov.fct <- function(std0, Cov0) sum(diag(std0%*%eval(Cov0)))
        trAsCov <- substitute(do.call(tr.fct, args=list(std0=std1, Cov0=Cov1)),
                              list(tr.fct = trAsCov.fct, std1=std, Cov1=Cov))
        r <- neighbor@radius
        Risk <- list(trAsCov = list(value = trAsCov,
                                    normtype = normtype),
                     asCov = Cov,
                     asBias = list(value = b, biastype = biastype,
                                   normtype = normtype,
                                   neighbortype = class(neighbor)),
                     asMSE = list(value = asMSE.0,
                                  r = r,
                                  at = neighbor))

        if(.withEvalAsVar) Risk <- .evalListRec(Risk)

        if(verbose)
           .checkPIC(L2deriv = L2deriv, neighbor = neighbor,
                     Distr = Distr, trafo = trafo, z = z, A = A, w = w,
                     z.comp = z.comp, A.comp = A.comp, ...)

        return(list(A = A, a = a, b = b, d = NULL, risk = Risk, info = info,
                    w = w, biastype = biastype, normtype = normtype,
                    call = mc, iter = iter, prec = prec, OIcall = OptIterCall))
    })




### helper function to check whether given b is in (bmin, bmax)
###        if not returns corresponding upper / lower case solution
.checkUpLow <- function(L2deriv, b, risk, neighbor, biastype, normtype,
                        Distr, Finfo, DistrSymm, L2derivSymm,
                        L2derivDistrSymm, z.start, A.start, trafo, maxiter,
                        tol, QuadForm, verbose, nrvalpts, warn){

            ClassIC <- trafo %*% solve(Finfo) %*% L2deriv

            lower.x <- getLow(Distr)
            upper.x <- getUp(Distr)
            x <- seq(from = lower.x, to = upper.x, length = nrvalpts)
            bmax <- sapply(x,function(x) evalRandVar(ClassIC,x))
            bmax <- sqrt(max(colSums(as.matrix(bmax^2))))
            cat("numerical approximation of maximal bound:\t", bmax, "\n")

            if(b >= bmax){
                if(warn) cat("'b >= maximum asymptotic bias' => (classical) optimal IC\n",
                             "in sense of Cramer-Rao bound is returned\n")
                res <- getInfRobIC(L2deriv = L2deriv, risk = asCov(), neighbor = neighbor,
                                    Distr = Distr, Finfo = Finfo, trafo = trafo,
                                    QuadForm = QuadForm, verbose = verbose)
                res <- c(res, list(biastype = biastype, normtype = normtype))
                trAsCov <- sum(diag(QuadForm%*%res$risk$asCov));
                r <- neighbor@radius
                res$risk$trAsCov <- list(value = trAsCov, normtype = normtype)
                res$risk$asBias <- list(value = b, biastype = biastype,
                                       normtype = normtype,
                                       neighbortype = class(neighbor))
                res$risk$asMSE <- list(value = trAsCov + r^2*b^2,
                                       r = r,
                                       at = neighbor)
                return(list(up=TRUE,low=FALSE,res=res))
            }

            res <- getInfRobIC(L2deriv = L2deriv,
                         risk = asBias(biastype = biastype, normtype = normtype),
                         neighbor = neighbor, Distr = Distr, DistrSymm = DistrSymm,
                         L2derivSymm = L2derivSymm, L2derivDistrSymm = L2derivDistrSymm,
                         z.start = z.start, A.start = A.start, trafo = trafo,
                         maxiter = maxiter, tol = tol, warn = warn, Finfo = Finfo,
                         verbose = verbose)
            bmin <- res$b

            cat("minimal bound:\t", bmin, "\n")
            if(b <= bmin){
                if(warn) cat("'b <= minimum asymptotic bias'\n",
                             "=> the minimum asymptotic bias (lower case) solution is returned\n")

                asMSE <- sum(diag(QuadForm%*%res$risk$asCov)) + neighbor@radius^2*bmin^2
                if(!is.null(res$risk$asMSE)) res$risk$asMSE <- asMSE
                   else     res$risk <- c(list(asMSE = asMSE), res$risk)

                return(list(up=FALSE,low=TRUE,res=res))
            }

            return(list(up=FALSE,low=FALSE,res=NULL))
}


################################################################################
## old routine
################################################################################
#
#setMethod("getInfRobIC", signature(L2deriv = "RealRandVariable",
#                                   risk = "asHampel",
#                                   neighbor = "UncondNeighborhood"),
#    function(L2deriv, risk, neighbor, Distr, DistrSymm, L2derivSymm,
#             L2derivDistrSymm, Finfo, trafo, onesetLM = FALSE,
#             z.start, A.start, upper = NULL, maxiter, tol, warn, verbose = FALSE,
#             checkBounds = TRUE){
#
#        biastype <- biastype(risk)
#        normtype <- normtype(risk)
#        p <- nrow(trafo)
#        if(! is(neighbor,"ContNeighborhood") && p>1)
#           stop("Not yet implemented")
#
#        FI <- solve(trafo%*%solve(Finfo)%*%t(trafo))
#        if(is(normtype,"InfoNorm") || is(normtype,"SelfNorm") )
#           {QuadForm(normtype) <- PosSemDefSymmMatrix(FI); normtype(risk) <- normtype}
#
#        std <- if(is(normtype,"QFNorm")) QuadForm(normtype) else diag(p)
#
#        if(is.null(z.start)) z.start <- numeric(ncol(trafo))
#        if(is.null(A.start)) A.start <- trafo
#        b <- risk@bound
#
#        if(checkBounds){
#            ClassIC <- trafo %*% solve(Finfo) %*% L2deriv
#            lower.x <- getLow(Distr)
#            upper.x <- getUp(Distr)
#            x <- seq(from = lower.x, to = upper.x, length = 5000)
#            bmax <- sapply(x,function(x) evalRandVar(ClassIC,x))
#            bmax <- sqrt(max(colSums(bmax^2)))
#            cat("numerical approximation of maximal bound:\t", bmax, "\n")
#            if(b >= bmax){
#                if(warn) cat("'b >= maximum asymptotic bias' => (classical) optimal IC\n",
#                             "in sense of Cramer-Rao bound is returned\n")
#                res <- getInfRobIC(L2deriv = L2deriv, risk = asCov(), neighbor = neighbor,
#                                    Distr = Distr, Finfo = Finfo, trafo = trafo,
#                                    QuadForm = std, verbose = verbose)
#                res <- c(res, list(biastype = biastype, normtype = normtype))
#                trAsCov <- sum(diag(std%*%res$risk$asCov));
#                r <- neighbor@radius
#                res$risk$trAsCov <- list(value = trAsCov, normtype = normtype)
#                res$risk$asBias <- list(value = b, biastype = biastype,
#                                       normtype = normtype,
#                                       neighbortype = class(neighbor))
#                res$risk$asMSE <- list(value = trAsCov + r^2*b^2,
#                                       r = r,
#                                       at = neighbor)
#                return(res)
#            }
#
#            res <- getInfRobIC(L2deriv = L2deriv,
#                         risk = asBias(biastype = biastype, normtype = normtype),
#                         neighbor = neighbor, Distr = Distr, DistrSymm = DistrSymm,
#                         L2derivSymm = L2derivSymm, L2derivDistrSymm = L2derivDistrSymm,
#                         z.start = z.start, A.start = A.start, trafo = trafo,
#                         maxiter = maxiter, tol = tol, warn = warn, Finfo = Finfo,
#                         verbose = verbose)
#            bmin <- res$b
#
#            cat("minimal bound:\t", bmin, "\n")
#            if(b <= bmin){
#                if(warn) cat("'b <= minimum asymptotic bias'\n",
#                             "=> the minimum asymptotic bias (lower case) solution is returned\n")
#
#                asMSE <- sum(diag(std%*%res$risk$asCov)) + neighbor@radius^2*bmin^2
#                if(!is.null(res$risk$asMSE)) res$risk$asMSE <- asMSE
#                   else     res$risk <- c(list(asMSE = asMSE), res$risk)
#
#                return(res)
#            }
#        }
#
#        comp <- .getComp(L2deriv, DistrSymm, L2derivSymm, L2derivDistrSymm)
#
#        z.comp <- comp$"z.comp"
#        A.comp <- comp$"A.comp"
#
#        z <- z.start
#        A <- A.start
#        if(is(neighbor,"ContNeighborhood")){
#            w <- new("HampelWeight")
#        }else if(is(neighbor,"TotalVarNeighborhood")){
#            w <- new("BdStWeight")
#        }
#        iter <- 0
#        a <- 0
#        repeat{
#            iter <- iter + 1
#            z.old <- z
#            A.old <- A
#
#            if(is(neighbor,"ContNeighborhood")){
#                clip(w) <- b
#                cent(w) <- z
#                stand(w) <- A
#            }else if(is(neighbor,"TotalVarNeighborhood")){
#                clip(w) <- if(iter==1) c(-b,b)/2 else c(0,b)+a
#                stand(w) <- A
#            }
#
#            weight(w) <- getweight(w, neighbor = neighbor, biastype = biastype,
#                                   normW = normtype)
#
#
#            z <- getInfCent(L2deriv = L2deriv, neighbor = neighbor,
#                            biastype = biastype, Distr = Distr, z.comp = z.comp,
#                            w = w)
#
#            if(is(neighbor,"TotalVarNeighborhood")){
#                  a <- z
#                  z <- solve(A,a)
#                  zc <- numeric(ncol(trafo))
#            }else if(is(neighbor,"ContNeighborhood")) {
#                  zc <- z
#            }
#
#            A <- getInfStand(L2deriv = L2deriv, neighbor = neighbor,
#                         biastype = biastype, Distr = Distr, A.comp = A.comp,
#                         cent = zc, trafo = trafo, w = w)
#
#            normtype.old <- normtype
#            if(is(normtype,"SelfNorm")){
#               normtype(risk) <- normtype <- updateNorm(normtype = normtype,
#                   L2 = L2deriv, neighbor = neighbor, biastype = biastype,
#                   Distr = Distr, V.comp = A.comp, cent = z, stand = A, w = w)
#               std <- QuadForm(normtype)
#            }
#
#            prec <- max(max(abs(A-A.old)), max(abs(z-z.old)))
#            if(verbose)
#                cat("current precision in IC algo:\t", prec, "\n")
#            if(prec < tol) break
#            if(iter > maxiter){
#                cat("maximum iterations reached!\n", "achieved precision:\t", prec, "\n")
#                break
#            }
#        }
#        if (onesetLM){
#            if(is(neighbor,"ContNeighborhood"))
#               cent(w) <- z
#            if(is(neighbor,"TotalVarNeighborhood"))
#               clip(w) <- c(0,b)+a
#            stand(w) <- A
#
#            weight(w) <- getweight(w, neighbor = neighbor, biastype = biastype,
#                                   normW = normtype)
#        }
#        else normtype <- normtype.old
#
#        info <- paste("optimally robust IC for 'asHampel' with bound =", round(b,3))
#        a <- as.vector(A %*% z)
#        Cov <- getInfV(L2deriv = L2deriv, neighbor = neighbor,
#                       biastype = biastype, Distr = Distr,
#                       V.comp = A.comp, cent = a,
#                       stand = A, w = w)
#        #getAsRisk(risk = asCov(), L2deriv = L2deriv, neighbor = neighbor,
#        #          biastype = biastype, Distr = Distr, clip = b, cent = a,
#        #          stand = A)$asCov
#        trAsCov <- sum(diag(std%*%Cov)); r <- neighbor@radius
#        Risk <- list(trAsCov = list(value = trAsCov,
#                                    normtype = normtype),
#                     asCov = Cov,
#                     asBias = list(value = b, biastype = biastype,
#                                   normtype = normtype,
#                                   neighbortype = class(neighbor)),
#                     asMSE = list(value = trAsCov + r^2*b^2,
#                                  r = r,
#                                  at = neighbor))
#
#        return(list(A = A, a = a, b = b, d = NULL, risk = Risk, info = info,
#                    w = w, biastype = biastype, normtype = normtype))
#    })










































