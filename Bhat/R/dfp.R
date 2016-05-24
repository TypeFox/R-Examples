"dfp" <-
function (x, f, tol=1e-5, nfcn = 0) 
{
	    #     Function Minimization for R. 
	    #     This function is part of the Bhat exploration tool and is
	    #     based on Fletcher's "Switching Method" (Comp.J. 13,317 (1970)).

	    #     This program is free software; you can redistribute it and/or modify
	    #     it under the terms of the GNU General Public License as published by
	    #     the Free Software Foundation; either version 2 of the License, or
	    #     (at your option) any later version.
	    #     This program is distributed in the hope that it will be useful,
	    #     but WITHOUT ANY WARRANTY; without even the implied warranty of
	    #     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	    #     GNU General Public License for more details.

            #     E Georg Luebeck (gluebeck@fhcrc.org)

            xt.inf <- 16
	    slamin <- 0.2
	    slamax <- 3
	    tlamin <- 0.05
	    tlamax <- 6
	    iter <- 0
            status <- 1
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
	    if (npar <= 0) {
	    	    warning("no. of parameters < 1")
	    	    stop() 
	    }
	    g <- numeric(npar)
	    gs <- numeric(npar)
	    g2 <- numeric(npar)
	    xt <- numeric(npar)
	    xxs <- numeric(npar)
	    dirin <- numeric(npar)
	    delgam <- numeric(npar)
	    vg <- numeric(npar)
	    ####  first call 
	    cat(date(), "\n","\n")
	    xt <- ftrf(x$est, x$low, x$upp)
	    fmin <- f(x$est)
	    nfcn <- nfcn + 1
            cat('starting at','\n')
	    cat(format(nfcn), "  fmin: ", fmin, "   ", format(x$est), "\n")
	    isw2 <- 0
	    nretry <- 0
	    nfcnmx <- 10000
	    vtest <- 0.001
	    apsi <- 0.05
	    ####  or .01 for better convergence?
	    up <- 1
	    ####  memorize current no. of function calls
	    npfn <- nfcn
	    rho2 <- 10 * apsi
	    rostop <- tol * apsi
	    trace <- 1
	    fs <- fmin
	    cat("\n")
	    # cat("rostop: ", rostop, "\n", "apsi:   ", apsi, "\n", "vtest: ", vtest, "\n")
	    dfp.loop <- 0
	    # while()
	    while (dfp.loop < 10) {
	    	    #### 1, 2 
	    	    #### ?? cat("COVARIANCE MATRIX NOT POSITIVE-DEFINITE","\n")
	    	    #### define step sizes dirin
                    d <- dqstep(list(label=x$label,est=btrf(xt, x$low, x$upp),low=x$low,upp=x$upp),f,sens=.01)
	    	    if (isw2 >= 1)  d <- 0.02 * sqrt(abs(diag(v)) * up)
	    	    dirin <- d

	    	    ##### obtain gradients, second derivatives, search for pos. curvature #########
	    	    ntry <- 0
	    	    negg2 <- 1
	    	    # loop 10

	    	    while (negg2 >= 1 & ntry <= 6) {
	    	    	    negg2 <- 0
	    	    	    #for(id in 1:npar)  loop 10
	    	    	    for (i in 1:npar) {
	    	    	    	    d <- dirin[i]
	    	    	    	    xtf <- xt[i]
	    	    	    	    xt[i] <- xtf + d
	    	    	    	    fs1 <- f(btrf(xt, x$low, x$upp))
	    	    	    	    nfcn <- nfcn + 1
	    	    	    	    xt[i] <- xtf - d
	    	    	    	    fs2 <- f(btrf(xt, x$low, x$upp))
	    	    	    	    nfcn <- nfcn + 1
	    	    	    	    xt[i] <- xtf
	    	    	    	    gs[i] <- (fs1 - fs2)/(2 * d)
	    	    	    	    g2[i] <- (fs1 + fs2 - 2 * fmin)/d^2
	    	    	    	    #if (g2[i]  <=  0.) 
	    	    	    	    if (g2[i] <= 0) {
	    	    	    	     ####  search if g2  <=  0. . .
	    	    	    	     cat("covariance matrix is not positive-definite or nearly singular", 
	    	    	    	      "\n")
	    	    	    	     negg2 <- negg2 + 1
	    	    	    	     ntry <- ntry + 1
	    	    	    	     d <- 50 * abs(dirin[i])
	    	    	    	     xbeg <- xtf
	    	    	    	     if (gs[i] < 0) {
	    	    	    	      dirin[i] <- -dirin[i]
	    	    	    	     }
	    	    	    	     kg <- 0
	    	    	    	     nf <- 0
	    	    	    	     ns <- 0
	    	    	    	     # while(ns < 10 ...)
	    	    	    	     while (ns < 10 & nf < 10) {
	    	    	    	      xt[i] <- xtf + d
	    	    	    	      f0 <- f(btrf(xt, x$low, x$upp))
	    	    	    	      nfcn <- nfcn + 1
	    	    	    	      #       cat("dfp search intermediate output:","\n")
	    	    	    	      #       cat("f0: ",f0,"  fmin: ",fmin,"   nfcn: ",nfcn,"\n")
	    	    	    	      #       cat("xt: ",xt,"\n")
	    	    	    	      if (xt[i] > xt.inf) {
                                        cat("parameter ", x$label[i], " near upper boundary","\n")}
	    	    	    	      if (xt[i] < -xt.inf) { 
	    	    	    	        cat("parameter ", x$label[i], " near lower boundary","\n")}
	    	    	    	      if (f0 == "NaN") {
	    	    	    	       warning("f0 is NaN")
	    	    	    	       nf <- 10
	    	    	    	       break
	    	    	    	      }
	    	    	    	      if (f0 <= fmin & abs(xt[i]) < xt.inf) {
	    	    	    	       # success	
	    	    	    	       xtf <- xt[i]
	    	    	    	       d <- 3 * d
	    	    	    	       fmin <- f0
	    	    	    	       kg <- 1
	    	    	    	       ns <- ns + 1
	    	    	    	      }
	    	    	    	      else {
	    	    	    	       if (kg == 1) {
	    	    	    	        ns <- 0
	    	    	    	        nf <- 0
	    	    	    	        # failure	
	    	    	    	        break
	    	    	    	       }
	    	    	    	       else {
	    	    	    	        kg <- -1
	    	    	    	        nf <- nf + 1
	    	    	    	        d <- (-0.4) * d
	    	    	    	       }
	    	    	    	      }
	    	    	    	     }
	    	    	    	     if (nf == 10) {
	    	    	    	      d <- 1000 * d
	    	    	    	      xtf <- xbeg
	    	    	    	      g2[i] <- 1
	    	    	    	      negg2 <- negg2 - 1
	    	    	    	     }
	    	    	    	     if (ns == 10) {
	    	    	    	      if (fmin >= fs) {
	    	    	    	       d <- 0.001 * d
	    	    	    	       xtf <- xbeg
	    	    	    	       g2[i] <- 1
	    	    	    	       negg2 <- negg2 - 1
	    	    	    	      }
	    	    	    	     }
	    	    	    	     xt[i] <- xtf
	    	    	    	     dirin[i] <- 0.1 * d
	    	    	    	     fs <- fmin
	    	    	    	    }
	    	    	    }
	    	    }
	    	    ####  provide output
	    	    if (ntry > 6) {
	    	    	    warning("could not find pos. def. covariance - procede with DFP")}
	    	    ntry <- 0
	    	    matgd <- 1
	    	    ####  diagonal matrix
	    	    ####  get sigma and set up loop
	    	    if (isw2 <= 1) {
	    	    	    ntry <- 1
	    	    	    matgd <- 0
	    	    	    v <- matrix(0, npar, npar)
	    	    	    diag(v) <- 2/g2
	    	    }
	    	    if (!all(diag(v) > 0)) {
	    	    	    ntry <- 1
	    	    	    matgd <- 0
	    	    	    v <- matrix(0, npar, npar)
	    	    	    # check whether always g2 > 0 ?
	    	    	    diag(v) <- 2/g2
	    	    }
	    	    xxs <- xt
	    	    sigma <- as.vector(gs %*% (v %*% gs)) * 0.5
	    	    if (sigma < 0) {
	    	    	    cat("covariance matrix is not positive-definite", 
	    	    	    	    "\n")
	    	    	    if (ntry == 0) {
	    	    	    	    #try one more time (setting ntry=1)
	    	    	    	    ntry <- 1
	    	    	    	    matgd <- 0
	    	    	    	    v <- matrix(0, npar, npar)
	    	    	    	    diag(v) <- 2/g2
	    	    	    	    xxs <- xt
	    	    	    	    sigma <- as.vector(gs %*% (v %*% gs)) * 0.5
	    	    	    	    ####  provide output
	    	    	    	    if (sigma < 0) {
	    	    	    	     isw2 <- 0
	    	    	    	     warning("could not find pos. def. covariance")
	    	    	    	     x$est <- btrf(xt, x$low, x$upp)
	    	    	    	     return(list(fmin = fmin, x = x))
	    	    	    	    }
	    	    	    }
	    	    }
	    	    else {
	    	    	    isw2 <- 1
	    	    	    iter <- 0
	    	    }
	    	    x$est <- btrf(xt, x$low, x$upp)
	    	    cat("dfp search intermediate output:", "\n")
	    	    cat("iter: ", iter, "  fmin: ", fmin, "   nfcn: ", 
	    	    	    nfcn, "\n")
	    	    ####  start main loop #######################################
	    	    if(dfp.loop == 0) {cat("start main loop:", "\n")} else {
                      cat("restart main loop:", "\n")}
	    	    f.main <- 0
	    	    #exit main if: f.main > 0
	    	    while (f.main == 0) {
	    	    	    #vector
	    	    	    dirin <- -0.5 * as.vector(v %*% gs)
	    	    	    #scalar
	    	    	    gdel <- as.vector(dirin %*% gs)
	    	    	    ####  linear search along -vg  . . .
	    	    	    xt <- xxs + dirin
	    	    	    xt[xt > xt.inf] <- xt.inf
	    	    	    xt[xt < -xt.inf] <- -xt.inf
	    	    	    f0 <- f(btrf(xt, x$low, x$upp))
	    	    	    nfcn <- nfcn + 1
	    	    	    ####  cat(format(nfcn),"   f0: ",f0,"   ",format(xt),"\n","\n")
	    	    	    ####  change to output on orig. scale
	    	    	    ####  quadr interp using slope gdel
	    	    	    denom <- 2 * (f0 - fmin - gdel)
	    	    	    if (denom <= 0) {
	    	    	    	    slam <- slamax
	    	    	    }
	    	    	    else {
	    	    	    	    slam <- -gdel/denom
	    	    	    	    if (slam > slamax) {
	    	    	    	     ####  cat("slam: ",format(slam),"\n")
	    	    	    	     slam <- slamax
	    	    	    	    }
	    	    	    	    else {
	    	    	    	     if (slam < slamin) 
	    	    	    	      slam <- slamin
	    	    	    	    }
	    	    	    }
	    	    	    # if (abs(slam-1.)  >=  0.1)
	    	    	    if (abs(slam - 1) >= 0.1) {
	    	    	    	    # else go to 70
	    	    	    	    xt <- xxs + slam * dirin
	    	    	    	    xt[xt > xt.inf] <- xt.inf
	    	    	    	    xt[xt < -xt.inf] <- -xt.inf
	    	    	    	    f2 <- f(btrf(xt, x$low, x$upp))
	    	    	    	    nfcn <- nfcn + 1
	    	    	    	    ####   cat(format(nfcn),"  f2: ",f2,"   ",format(xt),"\n")
	    	    	    	    ####   quadr interp using 3 points
	    	    	    	    aa <- fs/slam
	    	    	    	    bb <- f0/(1 - slam)
	    	    	    	    cc <- f2/(slam * (slam - 1))
	    	    	    	    denom <- 2 * (aa + bb + cc)
	    	    	    	    if (denom <= 0) {
	    	    	    	     tlam <- tlamax
	    	    	    	    }
	    	    	    	    else {
	    	    	    	     tlam <- (aa * (slam + 1) + bb * slam + cc)/denom
	    	    	    	     ####  cat("tlam: ",format(tlam),"\n"); cat(f0,f2,fs,'\n')
	    	    	    	     if (tlam > tlamax) {
	    	    	    	      tlam <- tlamax
	    	    	    	     }
	    	    	    	     else {
	    	    	    	      if (tlam < tlamin) 
	    	    	    	       tlam <- tlamin
	    	    	    	     }
	    	    	    	    }
	    	    	    	    xt <- xxs + tlam * dirin
	    	    	    	    xt[xt > xt.inf] <- xt.inf
	    	    	    	    xt[xt < -xt.inf] <- -xt.inf
	    	    	    	    f3 <- f(btrf(xt, x$low, x$upp))
	    	    	    	    nfcn <- nfcn + 1
	    	    	    	    if (f0 >= fmin & f2 >= fmin & f3 >= fmin) {
	    	    	    	     f.main <- 200
	    	    	    	     cat("break 200", "\n")
	    	    	    	     ####	cat(format(nfcn),"  f3: ",f3,"   ",format(xt),"\n")
	    	    	    	     break
	    	    	    	    }
	    	    	    	    if (f0 < f2 & f0 < f3) {
	    	    	    	     slam <- 1
	    	    	    	    }
	    	    	    	    else {
	    	    	    	     if (f2 < f3) {
	    	    	    	      f0 <- f2
	    	    	    	     }
	    	    	    	     else {
	    	    	    	      ## 55?	
	    	    	    	      f0 <- f3
	    	    	    	      slam <- tlam
	    	    	    	     }
	    	    	    	    }
	    	    	    	    dirin <- dirin * slam
	    	    	    	    xt <- xxs + dirin
	    	    	    	    xt[xt > xt.inf] <- xt.inf
	    	    	    	    xt[xt < -xt.inf] <- -xt.inf
	    	    	    }
	    	    	    fmin <- f0
	    	    	    isw2 <- 2
	    	    	    if (sigma + fs - fmin < rostop) {
	    	    	    	    f.main <- 170
	    	    	    	    ####  stop/convergence criteria
	    	    	    	    break
	    	    	    }
	    	    	    apsi.c <- sigma + rho2 + fs - fmin
	    	    	    if (apsi.c <= apsi) {
	    	    	    	    if (trace < vtest) {
	    	    	    	     f.main <- 170
	    	    	    	     break
	    	    	    	    }
	    	    	    }
	    	    	    #non-convergence
	    	    	    if (nfcn - npfn >= nfcnmx) {
	    	    	    	    f.main <- 190
	    	    	    	    break
	    	    	    }
	    	    	    iter <- iter + 1
	    	    	    ####  get gradient and sigma
	    	    	    ####  compute first and second (diagonal) derivatives 
	    	    	    fmin <- f(btrf(xt, x$low, x$upp))
	    	    	    nfcn <- nfcn + 1
	    	    	    ####      cat(format(nfcn),"  ",fmin,"   ",format(xt),"\n")
	    	    	    ####      cat("___________________________________________","\n")
	    	    	    for (i in 1:npar) {
	    	    	    	    vii <- v[i, i]
	    	    	    	    if (vii > 0) {
	    	    	    	     d <- 0.002 * sqrt(abs(vii) * up)
	    	    	    	    }
	    	    	    	    else {
	    	    	    	     d <- 0.02
	    	    	    	    }
	    	    	    	    xtf <- xt[i]
	    	    	    	    xt[i] <- xtf + d
	    	    	    	    fs1 <- f(btrf(xt, x$low, x$upp))
	    	    	    	    nfcn <- nfcn + 1
	    	    	    	    xt[i] <- xtf - d
	    	    	    	    fs2 <- f(btrf(xt, x$low, x$upp))
	    	    	    	    nfcn <- nfcn + 1
	    	    	    	    xt[i] <- xtf
	    	    	    	    g[i] <- (fs1 - fs2)/(2 * d)
	    	    	    	    g2[i] <- (fs1 + fs2 - 2 * fmin)/d^2
	    	    	    }
	    	    	    rho2 <- sigma
	    	    	    vg <- 0.5 * as.vector(v %*% (g - gs))
	    	    	    gvg <- as.vector((g - gs) %*% vg)
	    	    	    delgam <- as.vector(dirin %*% (g - gs))
	    	    	    sigma <- 0.5 * as.vector(g %*% (v %*% g))
	    	    	    if (sigma >= 0) {
	    	    	    	    if (gvg <= 0 | delgam <= 0) {
	    	    	    	     {
	    	    	    	      if (sigma < 0.1 * rostop) {
	    	    	    	       f.main <- 170
	    	    	    	       break
	    	    	    	      }
	    	    	    	      else {
	    	    	    	       f.main <- 1
	    	    	    	       break
	    	    	    	      }
	    	    	    	     }
	    	    	    	    }
	    	    	    }
	    	    	    else {
	    	    	    	    f.main <- 1
	    	    	    	    break
	    	    	    }
	    	    	    ####  update covariance matrix
	    	    	    trace <- 0
	    	    	    vii <- diag(v)
	    	    	    for (i in 1:npar) {
	    	    	    	    for (j in 1:npar) {
	    	    	    	     d <- dirin[i] * dirin[j]/delgam - vg[i] * 
	    	    	    	      vg[j]/gvg
	    	    	    	     v[i, j] <- v[i, j] + 2 * d
	    	    	    	    }
	    	    	    }
	    	    	    if (delgam > gvg) {
	    	    	    	    flnu <- dirin/delgam - vg/gvg
	    	    	    	    for (i in 1:npar) {
	    	    	    	     for (j in 1:npar) {
	    	    	    	      v[i, j] <- v[i, j] + 2 * gvg * flnu[i] * 
	    	    	    	       flnu[j]
	    	    	    	     }
	    	    	    	    }
	    	    	    }
	    	    	    xxs <- xt
	    	    	    gs <- g
	    	    	    trace <- sum(((diag(v) - vii)/(diag(v) + vii))^2)
	    	    	    trace <- sqrt(trace/npar)
	    	    	    fs <- f0
                          } # end main loop ###########################################
	    	    if (f.main == 1) {
	    	    	    x$est <- btrf(xt, x$low, x$upp)
	    	    	    dfp.loop <- dfp.loop + 1
	    	    	    cat("dfp loop: ", dfp.loop, "\n")
	    	    }
	    	    if (f.main == 190) {
	    	    	    #exceeds nfcn
	    	    	    isw1 <- 1
	    	    	    f.main <- 230
                            status <- 1
	    	    	    break
	    	    }
	    	    if (f.main == 200) {
	    	    	    cat("dfp search fails to converge, restart ...", 
	    	    	    	    "\n")
	    	    	    xt <- xxs
	    	    	    x$est <- btrf(xt, x$low, x$upp)
	    	    	    isw2 <- 1
	    	    	    cat("sigma: ", sigma, "\n")
	    	    	    if (sigma < rostop) {
	    	    	    	    f.main <- 170
	    	    	    }
	    	    	    else {
	    	    	    	    if (matgd >= 0) {
	    	    	    	     ###################################
	    	    	    	     dfp.loop <- dfp.loop + 1
	    	    	    	     cat("dfp loop: ", dfp.loop, "\n")
	    	    	    	    }
	    	    	    	    else {
                                      status <- 1
	    	    	    	      break
	    	    	    	    }
	    	    	    }
                     }
	    	    ####  CONVERGENCE
	    	    if (f.main == 170) {
                      status <- 0
	    	    	    cat("DFP search completed with status ",trunc(status),"\n","\n")
	    	    	    isw2 <- 3
	    	    	    if (matgd == 0) {
	    	    	    	    npargd <- npar * (npar + 5)/2
	    	    	    	    if (nfcn - npfn < npargd) {
	    	    	    cat("covariance matrix inaccurate - calculate Hessian","\n")
	    	    	    	     if (isw2 >= 2) {
	    	    	    	      isw2 <- 3
	    	    	    	      cat("perhaps dfp started near minimum - try newton-raphson", "\n")
                                    }
                          }
                                  }
	    	    	    break
	    	    }
                  }
	    fmin <- f(btrf(xt, x$low, x$upp)); nfcn <- nfcn + 1

	    x$est <- btrf(xt, x$low, x$upp)
            ####  compute error (logit scale)
            del <- dqstep(x,f,sens=.01)
            h <- logit.hessian(x,f,del,dapprox=FALSE,nfcn); nfcn <- h$nfcn
            v <- solve(h$ddf)
            
            xtl <- xt-1.96*sqrt(diag(v))
            xtu <- xt+1.96*sqrt(diag(v))
            
            ####  dfp search final output:
            xl <- btrf(xtl, x$low, x$upp)
            xu <- btrf(xtu, x$low, x$upp)
            
	    cat("Optimization Result:", "\n")
	    cat("iter: ", iter, "  fmin: ", fmin, "   nfcn: ", nfcn, 
	    	    "\n")
	    cat("\n")
	    m.out <- cbind(x$label, signif(x$est,8), signif(xl,8), signif(xu,8))
	    dimnames(m.out) <- list(1:npar, c("label", "estimate", "low", "high"))
	    print(m.out, quote = FALSE)
	    cat("\n")
	    return(list(fmin = fmin, label = x$label, est = x$est, status=status, nfcn=nfcn))
}




