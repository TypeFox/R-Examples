###################################################################################
# Lagrange Multipliers either by iteration or by optimization  --- new 10-08-09
###################################################################################

getLagrangeMultByIter <- function(b, L2deriv, risk, trafo,
                      neighbor, biastype, normtype, Distr,
                      a.start, z.start, A.start, w.start, std,
                      z.comp, A.comp, maxiter, tol,
                      verbose = NULL, warnit = TRUE){
        if(missing(verbose)|| is.null(verbose))
           verbose <- getRobAStBaseOption("all.verbose")
        LMcall <- match.call()

        ## initialization
        z <- z.start
        A <- A.start
        w <- w.start
        iter <- 0
        a <- a.start

        ## iteration-loop
        repeat{
            ## increment
            iter <- iter + 1
            a.old <- a
            z.old <- z
            A.old <- A

            ## update weight
            if(is(neighbor,"ContNeighborhood")){
                clip(w) <- b
                cent(w) <- as.numeric(z)
                stand(w) <- A
            }else if(is(neighbor,"TotalVarNeighborhood")){
                clip(w) <- if(.isVirginW(w)) c(-b,b)/2 else c(0,b)+a
                stand(w) <- A
            }
            weight(w) <- getweight(w, neighbor = neighbor, biastype = biastype,
                                   normW = normtype)

        #     print(w)
            ## update centering
            z <- getInfCent(L2deriv = L2deriv, neighbor = neighbor,
                            biastype = biastype, Distr = Distr, z.comp = z.comp,
                            w = w, tol.z = .Machine$double.eps^.5)
        #     print(c("z"=z))
            if(is(neighbor,"TotalVarNeighborhood")){
                  a <- z
                  z <- as.numeric(solve(A,a))
                  zc <- numeric(ncol(trafo))
            }else if(is(neighbor,"ContNeighborhood")) {
                  zc <- z
            }

            # update standardization
            A <- getInfStand(L2deriv = L2deriv, neighbor = neighbor,
                         biastype = biastype, Distr = Distr, A.comp = A.comp,
                         cent = zc, trafo = trafo, w = w)

        #     print(c("A"=A))
            ## in case of self-standardization: update norm
            normtype.old <- normtype
            if(is(normtype,"SelfNorm")){
               normtype(risk) <- normtype <- updateNorm(normtype = normtype,
                   L2 = L2deriv, neighbor = neighbor, biastype = biastype,
                   Distr = Distr, V.comp = A.comp, cent = zc, stand = A, w = w)
            }

            ## precision and iteration counting
            prec <- max(max(abs(A-A.old)), max(abs(a-a.old)),max(abs(z-z.old)))
#            if(verbose)
#              .checkPIC(L2deriv = L2deriv, neighbor = neighbor,
#                     Distr = Distr, trafo = trafo, z = zc, A = A, w = w,
#                     z.comp = z.comp, A.comp = A.comp)

            if(verbose && iter>1 && iter < maxiter && iter%%5 == 1){
                cat("current precision in IC algo:\t", prec, "\n")
                print(round(c(A=prec[1],a=prec[2]),3))
            }
            if(prec < tol) break
            if(iter > maxiter){
                if(warnit)
                   cat("maximum iterations reached!\n",
                       "achieved precision:\t", prec, "\n")
                break
            }
        }

        ## determine LM a
        if(is(neighbor,"ContNeighborhood"))
           a <- as.vector(A %*% zc)

        if(is(normtype,"QFNorm")) std <- QuadForm(normtype)

        return(list(A = A, a = a, z = zc, w = w,
                    biastype = biastype, normtype = normtype,
                    normtype.old = normtype.old,
                    risk = risk, std = std,
                    iter = iter, prec = prec, b = b,
                    call = LMcall ))
}

getLagrangeMultByOptim <- function(b, L2deriv, risk, FI, trafo,
                      neighbor, biastype, normtype, Distr,
                      a.start, z.start, A.start, w.start, std, z.comp,
                      A.comp, maxiter, tol, verbose = NULL, ...){

        if(missing(verbose)|| is.null(verbose))
           verbose <- getRobAStBaseOption("all.verbose")
        LMcall <- match.call()
        ### manipulate dots in call -> set control argument for optim
        dots <- list(...)
        if(is.null(dots$method)) dots$method <- "L-BFGS-B"

        if(!is.null(dots$control)){
            if(is.null(dots$control$maxit)) dots$control$maxit <-  round(maxiter)
            if(is.null(dots$control$reltol)) dots$control$reltol <- tol
            if(is.null(dots$control$abstol)) dots$control$abstol <- tol
        }else{
            dots$control = list(maxit=min(round(maxiter),1e8), reltol=tol, abstol=tol)
        }
        #print(dots$control)
        ## initialization
        z <- z.start
        A <- A.start
        p <- nrow(trafo)
        k <- ncol(trafo)

        A0vec0 <- as.numeric(cbind(A, A%*%z))

        lvec0 <- seq(along=A0vec0)
        A0log <- as.logical(cbind(A.comp, as.logical(A.comp%*%as.numeric(z.comp)>0)))
        lvlog <- lvec0[A0log]
        A0vec1 <- A0vec0[A0log]
#        print(list(A0vec0,A0log,lvlog,A0vec1))

        iter1 <- 0
        stdC  <- stdC.opt <- std
        optV <- Inf
        
        risk1 <- risk1.opt <- risk
        normtype1 <- normtype1.old <- normtype
        normtype1.opt <- normtype1.opt.old <- normtype
        w1.opt <- w1 <- w.start
        z1.opt <- numeric(k)
        b.opt <- b

        optimfct <- function(A0vec){
            iter1 <<- iter1 + 1
#            print(A0vec)
            A0vecA <- numeric(p*(k+1))

            A0vecA[lvlog] <- A0vec

            ### read out current value of LM in usual format
            A0 <- matrix(A0vecA[1:(p*k)],nrow=p,ncol=k)
            a0 <- as.numeric(A0vecA[(p*k)+(1:p)])

#            print(list(A0vecA,A0,a0))

            z0 <- as.numeric(solve(A0,a0))
            std0 <- stdC
            w0 <- w1
            risk0 <- risk1
            b0 <- b
            
            if(is(risk0,"asMSE")){
            funint.opt <-
                   function(b1){
                      -getInfGamma(L2deriv = L2deriv, risk = risk0,
                                 neighbor = neighbor, biastype = biastype,
                                 Distr = Distr, stand = A0, cent = z0, clip = b1,
                                 power = 2)+radius(neighbor)^2*b1^2
                      }

            b0 <- optimize(funint.opt, interval=c(1e-8,1e8))$minimum
            }

            ### determine corresponding weight
            if(is(neighbor,"ContNeighborhood")){
                clip(w0) <- b0
                cent(w0) <- as.numeric(z0)
                stand(w0) <- A0
            }else if(is(neighbor,"TotalVarNeighborhood")){
                clip(w0) <- if(.isVirginW(w0)) c(-b0,b0)/2 else c(0,b0)+a0
                stand(w0) <- A0
            }
            weight(w0) <- getweight(w0, neighbor = neighbor, biastype = biastype,
                                    normW = normtype1)

            ### in case of self-standardization update norm:
            if (is(normtype1,"SelfNorm"))
                {
                   ## transport current & precedent norm outside the optimizer:
                   normtype1.old <<- normtype1
                   normtype1 <<- updateNorm(normtype = normtype1,
                                            L2 = L2deriv, neighbor = neighbor,
                                            biastype = biastype, Distr = Distr,
                                            V.comp = A.comp, cent = a0,
                                            stand = A0, w = w0)
                   normtype(risk0) <- normtype1
                   ## transport current quadratic form & risk outside
                   ## the optimizer:
                   std0 <- stdC <<- QuadForm(normtype1)
                   risk1 <<- risk0
                   }

            ### for Y_A = A Lambda - a and D
            ### use LM A1 = (A,a) is minimizer of
            ### (*=c, arbitrary Biasterm):
            ##      A1 -> E |Y_{A1}|_Q^2 / 2 + gamma(Q,A1,b)/2 - tr QAD'
            ### (*=v, p=1):
            ##      A1 -> E Y_{A1}^2/2 + gamma(A1,b)/2 - tr AD' - E Y_{A1,-}^2/2
            ###### for gamma([Q,]A,b) = E[{Y_A (1-w_b(|Y_A|_Q))}^2]

            riskA <- risk0
            if(is(riskA, "asHampel")){
               riskA <- asMSE(biastype=biastype, normtype=normtype)
               val <-   (as.numeric(t(a0)%*%std0%*%a0)/2 +
                          sum(diag(std0%*%A0%*%FI%*%t(A0)))/2 +
                          ## ~ E |Y_A|_Q^2 / 2
                          getInfGamma(L2deriv = L2deriv, risk = riskA,
                                 neighbor = neighbor, biastype = biastype,
                                 Distr = Distr, stand = A0, cent = z0, clip = b0,
                                 power = 2)/2 -
                           # ~ - E[|Y_A|_Q^2 (1-w_b(|Y_A|_Q))^2]/2
                           sum(diag(std0%*%A0%*%t(trafo)) ))
                        ## ~tr_Q AD'

               ## in case TotalVarNeighborhood additional correction term:
               if(is(neighbor,"TotalVarNeighborhood"))
                  val <- (val -a0^2/2 -
                          E(Distr, fun = function(x){ ## ~ - E Y_-^2/2
                              L2 <- evalRandVar(L2deriv, as.matrix(x)) [,,1]- z0
                              Y <- A0 %*% L2
                              return(Y^2*(Y<0))
                              },  useApply = FALSE)/2)

            }else if(is(riskA,"asMSE")){
               val <- (E(object = Distr, fun = function(x){
                          X <- evalRandVar(L2deriv, as.matrix(x))[,,1] - z0
                          Y <- A0 %*% X
                          nY <- norm(risk0)(Y)
                          return(nY^2*weight(w0)(X))
                        },  # E|Y|^2 w
                        useApply=FALSE) /2 -
                       sum(diag(std0%*%A0%*%t(trafo)) ))
                     ## ~tr_Q AD'

               ## in case TotalVarNeighborhood additional correction term:
               if(is(neighbor,"TotalVarNeighborhood"))
                  val <- (val -a0^2/2 -
                          E(Distr, fun = function(x){
                              X <- evalRandVar(L2deriv, as.matrix(x))[,,1] - z0
                              Y <- A0 %*% X
                              return(Y^2*(Y<0))
                             },
                            useApply=FALSE)/2)
            }

            ## if this is the current optimum
            ## transport some values outside the optimizer:
#            print(val)
            if(val < optV){
#            print(list("INNEN",A0=A0,a0=a0,b=b,z0=z0,z01=MASS::ginv(A0)%*%a0,val=val,"ENDINNEN"))
               optV <<- val
               w1.opt <<- w0
               z1.opt <<- z0
               b.opt <<- b0
               normtype1.opt.old <<- normtype1.opt
               normtype1.opt <<- normtype1
               risk1.opt <<- risk0
               stdC.opt <<- stdC
            }

            if(verbose && iter1>1 && iter1 < maxiter && iter1%%8 == 1){
                print(round(c(val=val,A=A0,a=a0,b=b0),4))
            }

            return(val)
        }

        ### optimization:
        opterg <- do.call(optim,
                          args = c(list(par = A0vec1, fn = optimfct),dots))

#        print(opterg)
        ### read out results of call to "optim"
        Aoptvec <- numeric(p*(k+1))
        Aoptvec[lvlog] <- opterg$par

        A <- matrix(Aoptvec[1:(p*k)],nrow=p,ncol=k)
        a <- as.numeric(Aoptvec[(p*k)+(1:p)])
        z <- z1.opt
        w <- w1.opt
        b <- b.opt

        return(list(A = A, a = a, z = z, w = w,
                    biastype = biastype, normtype = normtype1.opt,
                    normtype.old = normtype1.opt.old,
                    risk = risk1.opt, std = stdC.opt, iter = iter1,
                    prec = opterg$convergence, b = b,
                    call = LMcall ))
}

