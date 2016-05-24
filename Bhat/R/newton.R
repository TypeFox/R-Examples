"newton" <-
function (x, f, eps=1e-1, itmax=10, relax=0, nfcn = 0) 
{
	    #     General Newton-Raphson module written in R. 
	    #     This module is part of the Bhat likelihood exploration tool.

	    #     This program is free software; you can redistribute it and/or modify
	    #     it under the terms of the GNU General Public License as published by
	    #     the Free Software Foundation; either version 2 of the License, or
	    #     (at your option) any later version.
	    #     This program is distributed in the hope that it will be useful,
	    #     but WITHOUT ANY WARRANTY; without even the implied warranty of
	    #     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	    #     GNU General Public License for more details.

            #     E. Georg Luebeck, Fred Hutchinson Cancer Research Center
	    #     Fairview Avenue N, Seattle, WA 98109-1024

	    if (!is.list(x)) {
	    	    cat("x is not a list! see help file", "\n")
	    	    return()
	    }
	    names(x)[1] <- "label"
	    names(x)[2] <- "est"
	    names(x)[3] <- "low"
	    names(x)[4] <- "upp"
	    npar <- length(x$est
                           )
	    ####  objects:
	    if (npar <= 0) {
	    	    warning("no. of parameters < 1")
	    	    stop() 
	    }

            small <- 1.e-8
            nlm <- (npar*npar+3*npar)/2
            np2 <- 2*npar
            sens <- 0.1
            flag.hessian <- 1

	    ####  first call 
	    cat(date(), "\n")
	    xt <- ftrf(x$est, x$low, x$upp)
	    f0 <- f(x$est); nfcn <- nfcn + 1

            for(n in 1:itmax) {

              cat(nfcn, '  fmin: ', format(f0), "   ", format(x$est), '\n')

                # *** INITIALIZE NEWTON BY CALCULATION OF A
                #     NLM POINT SIMPLEX GEOMETRY IN NPAR DIM PARAMETER SPACE
                xn <- matrix(rep(xt,nlm),ncol=nlm)
                # ***   ON AXES
                if(n==1) del <- dqstep(x,f,sens)
                if(flag.hessian==1) {
                  h <- logit.hessian(x,f,del,dapprox=FALSE,nfcn)  # returns list: $df,$ddf,$nfcn
                  nfcn <- h$nfcn
                  ddf <- h$ddf; df <- h$df
                } else {

                # ***	compute delta distance in each direction
                for(i in 1:npar) {
                  j.even <- 2*i
                  j.odd  <- j.even-1
                  xn[i,j.even] <- xn[i,j.even]-del[i]
                  xn[i,j.odd]  <- xn[i,j.odd] +del[i]
                }

                # ***   OFF AXIS
                if(npar > 1) {
                  mc <- np2+1
                  for(i in 2:npar) {
                    for(j in 1:(i-1)) {
                      xn[i,mc] <- xn[i,mc]+del[i]
                      xn[j,mc] <- xn[j,mc]+del[j]
                      mc <- mc+1
                    }}
                }

                f.vec <- numeric(nlm) 
                for(i in 1:nlm) {
                  f.vec[i] <- f(btrf(xn[,i], x$low, x$upp))
	    	  nfcn <- nfcn + 1
                }
                
                # ***   FIRST AND DIAGONAL SECOND DERIVATIVES
                i <- 1:npar
                i.even <- 2*i
                i.odd <- i.even-1

                df <- (f.vec[i.even]-f.vec[i.odd])/2/del
                ddf <- diag((f.vec[i.even]+f.vec[i.odd]-2*f0)/del/del/2)

                # ***   SECOND DERIVATIVES

                if(npar > 1) {
                  mc <- np2+1
                  for(i in 2:npar) {
                    for(j in 1:(i-1)) {
                      ip <- 2*i-1; jp <- 2*j-1
                      ddf[i,j] <- (f0+f.vec[mc]-f.vec[jp]-f.vec[ip])/del[i]/del[j]
                      mc <- mc+1
                    }}

                  ddf <- ddf+t(ddf)
                }
              
                # ***   EIGEN VALUES and MATRIX INVERSION

                eig <- eigen(ddf)
                cat('approx. eigen values: ',format(eig$values),'\n')
                if(any(eig$values < 0)) {
                  warning('approx. hessian not pos. definite')
                }
                if(any(abs(eig$values) < small)) {
                  warning('approx. hessian may be singular')
                }
              }
                # cat('compare:','\n')
                # print(format(h$ddf),quote=FALSE)
                # print(format(ddf),quote=FALSE)
                      
                b <- diag(1,npar)
                # *** if inversion fails need to return or reinitialize NR (not implemented)
                ggf <- solve(ddf,b,tol=1.e-10)

                # cat('unity test:','\n'); print(format(ggf %*% ddf),quote=FALSE)
              
                # ---	variable Newton stepping

                disc <- ggf  %*% df
                disc0 <- disc
                rneps <- 1; nepsc <- 0; xt0 <- xt; f1 <- f0

                xt <- xt0 + disc
                # cat(xt0,'\n',xt,'\n',df,'\n')
                tmp <- numeric(npar)
                tmp[abs(df) <= eps] <- 1
                nz <- sum(tmp)
                f0 <- f(btrf(xt, x$low, x$upp)); nfcn <- nfcn + 1

                if(relax > 0) {
                while(nepsc < 8 && f0 > f1 ) {
                  rneps <- .75*rneps
                  cat('reducing step size: ',format(rneps,f0),'\n')
                  nepsc <- nepsc+1
                  disc <- rneps * disc0
                  xt <- xt0 + disc
                  f0 <- f(btrf(xt, x$low, x$upp)); nfcn <- nfcn + 1
                 }
                }

                if(nepsc == 8) warning('problem finding lower function value, will continue')
                
                # ---   STOP-CRITERION FOR ITERATION --------------------------------------
                if(nz == npar && min(h$eigen) > small) {

                status <- 'converged'

                # ---	compute Hessian symmetrically:
                x$est <- btrf(xt,x$low,x$upp)  
                # h <- logit.hessian(x,f,del,dapprox=FALSE,nfcn)  # returns list: $df,$ddf,$nfcn
                b <- diag(1,npar)
                ggf <- solve(ddf,b,tol=1.e-10)

                xtl <- rep(NA,npar); xtu <- rep(NA,npar)
                se <- numeric(npar)
                ggf.diag <- diag(ggf)
                se[ggf.diag >= 0] <- sqrt(ggf.diag[ggf.diag >= 0])
                xtl <- xt - 1.96 * se
                xtu <- xt + 1.96 * se

                x$est <- btrf(xt, x$low, x$upp)
                xlow  <- btrf(xtl, x$low, x$upp)
                xupp  <- btrf(xtu, x$low, x$upp)

                cat('\n')  
                cat('Bhat run: ',date(),' status: ',status,'\n')
                # cat('Optimization Result:', '\n')
                cat("iter: ", n, "  fmin: ", format(f0), "   nfcn: ", nfcn, "\n")
                cat("\n")
                m.out <- cbind(x$label, round(x$est,6), round(h$df,6), round(xlow,6), round(xupp,6))
                dimnames(m.out) <- list(1:npar, c("label", "estimate", "log deriv", "lower 95%" , "upper 95%"))

                print(m.out, quote=FALSE)
                cat("\n")
                return(list(fmin = f0, label = x$label, est = x$est, low = xlow, upp = xupp))

              }

                # --- graphical diagnostic (do be done)

                # ***	DURING NON-CONVERGEING CYCLES
                status <- 'non-converged'
                xtl <- rep(NA,npar); xtu <- rep(NA,npar)
                se <- numeric(npar)
                ggf.diag <- diag(ggf)
                se[ggf.diag >= 0] <- sqrt(ggf.diag[ggf.diag >= 0])
                xtl <- xt - 1.96 * se
                xtu <- xt + 1.96 * se

                x$est <- btrf(xt, x$low, x$upp)
                xlow  <- btrf(xtl, x$low, x$upp)
                xupp  <- btrf(xtu, x$low, x$upp)

                if(n == itmax) {
                  cat('\n')
                  cat('Bhat run: ',date(),' status: ',status,'\n')
                  # cat('Optimization Result:', '\n')
                } 
              
                cat('\n',"iter: ", n, "  fmin: ", f0, "   nfcn: ", nfcn, 
          	    "\n")
                cat("\n")
                vbar <- rep("|",npar)
                m.out <- cbind(x$label, round(x$est,6), round(df,6), round(xlow,6), round(xupp,6))
                dimnames(m.out) <- list(1:npar, c("label", "estimate", "log deriv", "lower 95%" , "upper 95%"))
                print(m.out, quote = FALSE)
                cat("\n")
            }
                return(list(fmin = f0, label = x$label, est = x$est, low = xlow, upp = xupp))
        }
