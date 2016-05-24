"mymcmc" <-
function (x, nlogf, m1, m2=m1, m3, scl1=0.5, scl2=2, skip=1, covm=0, nfcn = 0, plot=FALSE)
{
	    #     MCMC/MH sampler for R
	    #     This module is part of the Bhat likelihood exploration tool.

	    #     This program is free software; you can redistribute it and/or modify
	    #     it under the terms of the GNU General Public License as published by
	    #     the Free Software Foundation; either version 2 of the License, or
	    #     (at your option) any later version.
	    #     This program is distributed in the hope that it will be useful,
	    #     but WITHOUT ANY WARRANTY; without even the implied warranty of
	    #     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	    #     GNU General Public License for more details.

            #     E. Georg Luebeck (gluebeck@fhcrc.org)

	    if (!is.list(x)) {
	    	    cat("x is not a list! see help file", "\n")
	    	    return()
	    }
	    names(x)[1] <- "label"
	    names(x)[2] <- "est"
	    names(x)[3] <- "low"
	    names(x)[4] <- "upp"
	    npar <- length(x$est)
            
	    ####  objects:
	    if (npar <= 0) stop('no. of parameters must be >= 1')

            small <- 1.e-8
            sens <- 0.1
            xinf <- 20.

            # initialize graphical output
            if(plot==TRUE) {
            par(mfrow=c(npar+1,1),bg="grey")}

            # *** initialize f.old
	    cat(date(), "\n")
	    xt <- ftrf(x$est, x$low, x$upp)
	    f.old <- nlogf(x$est); nfcn <- nfcn + 1
            f.first <- f.old

            acc.count <- 0

            # *** if COVM is not provided
            if(!is.matrix(covm)) {

            # initialize monitors
            f.mon <- numeric(m1+m3); x.mon <- matrix(0,ncol=npar,nrow=m1+m3)

            cat('trying a proposal COVM using the Hessian','\n')
            del <- dqstep(x,nlogf,sens)
            h <- logit.hessian(x,nlogf,del,dapprox=FALSE,nfcn)  # returns list: $df,$ddf,$nfcn
            nfcn <- h$nfcn
            ddf <- h$ddf

            # ***   EIGEN VALUES and MATRIX INVERSION
            eig <- eigen(ddf)
            
            if(any(eig$values < small)) {
                  cat('Hessian may not be positive definite','\n')
                  cat('trying a proposal COVM using dqstep','\n')
                  del <- dqstep(x,nlogf,4.)
                  eig <- .1/del/del; ddf <- diag(eig); eig <- eigen(ddf)
                } 

            # ggf <- 0.5 * solve(ddf,diag(1,npar),tol=1.e-10)
            
            # cat('unity test:','\n'); print(format(ggf %*% ddf),quote=FALSE)
            cat('first pilot chain: using MLEs and log-transformed covariance','\n','\n')
            
            nc <- 0
            for (n in 1:m1) {

              accept <- 1

              # *** compute proposal x' (=yt). If unacceptable, reject
              dx <-  eig$vectors %*% rnorm(npar,0,1/(0.5*sqrt(eig$values)))
              yt <- xt + scl1 * dx
                                          
              # *** get log-likelihood
              f.new <- nlogf(btrf(yt, x$low, x$upp)); nfcn <- nfcn + 1
                                      
              # *** boundary checks from within func ...

              # *** R ratio 
              if(any(abs(yt) > xinf)) {
                cat('cycle',n,': mymcmc close to boundary, move rejected!','\n'); accept <- 0} else {
              accept <- min(accept,exp(-f.new+f.old))}
              
              if(accept == 1) {xt <- yt; f.old <- f.new; acc.count <- acc.count+1
                             xt.maxl <- yt; f.maxl <- f.new # approximate/search max-likelihood
                             } else {
                               if(runif(1) <= accept) {
                                         xt <- yt; f.old <- f.new; acc.count <- acc.count+1}
                             }

              f.mon[n] <- f.old; x.mon[n,] <- btrf(xt, x$low, x$upp)

              # *** regular output and graphical monitor
              if(n%%(100*skip) == 0) {
                m.out <- c("n:",signif(n,3),"acceptance:",signif(acc.count/100/skip,3),"-log lkh:",signif(f.new,8),signif(x.mon[n,],6))
                cat(m.out,'\n')
                acc.count <- 0
                n.skip <- seq(skip,n,skip)

                if(plot==TRUE) {
                par(mar=c(0, 5, 0.1, 4) + 0.1)
                plot(f.mon[n.skip], type='l', xlab = " ", ylab = "-log-density",col=2)
                for (i in 1:(npar-1)) {
                  par(mar=c(0, 5, 0, 4) + 0.1)
                  plot(x.mon[n.skip,i], type='l', xlab = " ", ylab = x$label[i], col=2)
                }
                par(mar=c(0.1, 5, 0, 4) + 0.1)
                plot(x.mon[n.skip,npar], type='l', xlab = " ", xaxt='n', ylab = x$label[npar], col=2)
                nc <- nc + 1
              }
              }
            }

            # *** obtain empirical covariance of increments
            covm <- cov(x.mon[2:m1,]-x.mon[1:(m1-1),]) # now redundant
            cat('\n','\n')
            # print(covm,quote=FALSE)
            # cat('\n','\n')

            # *** pilot run 2 *** includes m2 cycles for incremental adjustment of covm
            cat('second pilot run:','\n','\n')
            eig <- eigen(covm)
            cat('eigen values:',eig$values,'\n','\n')

          } else {
            # covm for mvn proposal distribution given as input  
            if(nrow(covm)!=length(x$est)) {stop('dimension of covm not specified correctly')}
            nc <- 0 
            # re-initialize monitors
            m1 <- 1
            f.mon <- numeric(m1+m3); x.mon <- matrix(0,ncol=npar,nrow=m1+m3)
            eig <- eigen(covm); x.mon[1,] <- x$est; f.mon <- f.first
          }
            acc.count <- 0
            for (n in (m1+1):(m1+m3)) {

              x.mon[n,] <- x.mon[n-1,]; f.mon[n] <- f.mon[n-1]
              accept <- 1

              # *** compute proposal x' (=yt). If unacceptable, reject
              dx <-  eig$vectors %*% rnorm(npar,0,sqrt(eig$values))
              y <- x.mon[n-1,] + scl2 * dx
              if(any((y-x$low) < small)) {cat('move rejected (lower bound)','\n'); accept <- 0}
              if(any((x$upp-y) < small)) {cat('move rejected (upper bound)','\n'); accept <- 0}

              # *** boundary checks from within func ...

              if(accept > 0) {
              # *** get log-likelihood
              f.new <- nlogf(y); nfcn <- nfcn + 1
                                      
              # *** acceptance ratio and updates
              accept <- min(accept,exp(-f.new+f.mon[n]))
              if(accept == 1) {x.mon[n,] <- y; f.mon[n] <- f.new; acc.count <- acc.count+1
                             x.maxl <- y; f.maxl <- f.new # approximate/search max-likelihood
                             } else {
                               if(runif(1) <= accept) {
                                 x.mon[n,] <- y; f.mon[n] <- f.new; acc.count <- acc.count+1
                               }
                             }
            }

              # *** regular output and graphical monitor
              if(n%%(100*skip) == 0) {
                m.out <- c("n:",signif(n,3),"acceptance:",signif(acc.count/100/skip,3),"-log lkh:",signif(f.new,8),signif(x.mon[n,],6))
                cat(m.out,'\n')
                acc.count <- 0

                n.skip.1 <- seq(skip,n,skip)
                n.skip.2 <- seq(skip,min(n,m1+m2),skip)

                if(plot==TRUE) {
                if(m1 > 1) {brncol <- 3} else {brncol <- 2}
                  
                par(mar=c(0, 5, 0.1, 4) + 0.1)
                plot(f.mon[n.skip.1], type='l', xlab = " ", ylab = "-log-density",col=2)
                lines(f.mon[n.skip.2],col=brncol) #pilot cycles
                for (i in 1:(npar-1)) {
                  par(mar=c(0, 5, 0, 4) + 0.1)
                  plot(x.mon[n.skip.1,i], type='l', xlab = " ", ylab = x$label[i], col=2)
                  lines(x.mon[n.skip.2,i],col=brncol) #pilot cycles
                }
                par(mar=c(0.1, 5, 0, 4) + 0.1)
                plot(x.mon[n.skip.1,npar], type='l', xlab = " ", xaxt='n', ylab = x$label[npar], col=2)
                lines(x.mon[n.skip.2,npar],col=brncol) #pilot cycles
                
                nc <- nc + 1
              }

                if(m1 > 1) { #note: when covm is passed to mymcmc, m1 is set to 1
                # update covariance using sampled increments (m1+1):n
                if(n <= m1+m2) {
                covm <- cov(x.mon[2:n,]-x.mon[1:(n-1),])
                eig <- eigen(covm)
                if(any(eig$values < small)) {warning('covm nearly singular')}
              }
              }
              }
            }
            return(list(f=f.mon,mcmc=x.mon,covm=covm))
            }











