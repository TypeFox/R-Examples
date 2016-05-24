'arima' <-function (x, order = c(0, 0, 0), seasonal = list(order = c(0, 
    0, 0), period = NA), xreg = NULL, include.mean = TRUE, transform.pars = TRUE, 
    fixed = NULL, init = NULL, method = c("CSS-ML", "ML", "CSS"), 
    n.cond, optim.control = list(), kappa = 1e+06, io = NULL, 
    xtransf, transfer = NULL) 
{
#
#  This function is based on the arima function of the stats package
#  of R. Below the copright statement of the arima function is reproduced. 
#
#  File src/library/stats/R/arima.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 2002-12 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

#
# (arimax) Date: April, 18, 2006
# Revised: August 15, 2006; Oct, 26, 2012; Nov 12, 2012


deriv=function(u){
#
# compute the derivative of the transformation from u to phi given on p. 128
# of Jones (1993) "Longitudinal Data with Serial Correlation: a state-space approach.
#
p=length(u)
#if(p<=1) return(2*exp(-u)/(1+exp(-u))^2)
#A=diag(as.vector(2*exp(-u)/(1+exp(-u))^2))
#pacf=(1-exp(-u))/(1+exp(-u))

if(p<=1) return(4*exp(2*u)/(1+exp(2*u))^2)
A=diag(as.vector(4*exp(2*u)/(1+exp(2*u))^2))
pacf=tanh(u)

init=diag(as.vector(rep(1,p)))
fd=NULL
for ( i in 1:p) {
oldphi=pacf[1]
oldderiv=init[i,]
for (k in 2:p) {
mult=pacf[k]
newphi=rep(NA,k)
newderiv=rep(NA,k)
newphi[k]=pacf[k]
newderiv[k]=init[i,k]
mult1=newderiv[k]
for (j in 1: (k-1)) {
newphi[j]=oldphi[j]-mult*oldphi[k-j]
newderiv[j]=oldderiv[j]-mult*oldderiv[k-j]-mult1*oldphi[k-j]
}
oldphi=newphi
oldderiv=newderiv
}
fd=rbind(fd, newderiv)
}
A%*%fd
}

difpar=function(par,arma){
#
# compute the first derivative matrix of the transformation from (u,U) to 
# (phi, Phi).
#
res=diag(length(par))
 if(arma[1]) res[1:arma[1],1:arma[1]]=(deriv(par[1:arma[1]]))
 if(arma[3]) {indx=sum(arma[c(1,2)])+seq(arma[3])
              res[indx,indx]=(deriv(par[indx]))}
res
}

undotrans=function(par,arma){
#
# compute the transformation from (u,U) to (phi,Phi).
#
      if(arma[1]) par[1:arma[1]]=levinson(par[1:arma[1]])
      if(arma[3]) par[sum(arma[c(1,2)])+seq(arma[3])]=levinson(par[sum(arma[c(1,2)])+seq(arma[3])])
par
}

pacfparm=function(phi){
#
# compute the transformation from phi to u.
#
res=stats:::ARMAacf(ar=phi,lag.max=length(phi),pacf=TRUE)
#log((1+res)/(1-res))
atanh(res)
}

transformpar=function(par,arma){
#
# conpute the transformation from (phi,Phi) to (u,U).
#
 if(arma[1]) par[1:arma[1]]=pacfparm(par[1:arma[1]])
 if(arma[3]) par[sum(arma[c(1,2)])+seq(arma[3])]=pacfparm(par[sum(arma[c(1,2)])+seq(arma[3])])
 par
}
 

conv=function(x,y){m=length(x)+length(y)-1
         res=rep(0,m)
         fres=res
         n=length(x)
         indx=1:n
         for (i in y) {
            res1=res
            res1[indx]=x*i
            fres=fres+res1
            indx=indx+1
          }
fres
}
    conv1=function(x, period=1){c(x,rep(0, period))-c(rep(0, period),x)}
    if (!missing(transfer)) {
        for (i in seq(length(transfer))) transfer[[i]][2] = transfer[[i]][2] + 
            1
    }
    checkIO = function(x, phi = NULL, theta = NULL, Delta = NULL) {
        name1 = colnames(x)
        d = max(length(phi) + length(Delta), length(theta)) + 
            1
        ind = 0
        for (i in name1) {
            ind = ind + 1
            if (length(grep("IO-", i)) >= 1) 
                x[, ind] = armafilter(c(rep(0, d), x[, ind]), 
                  phi = phi, theta = c(1, theta), Delta = Delta)[-(1:d)]
        }
        x
    }
    armafilter = function(x, phi, theta, Delta) {
        if (length(Delta) > 0) 
            x = filter(x, filter = Delta, sides = 1, method = "recursive")
        if (length(phi) > 0 && any(phi != 0)) 
            x = filter(x, filter = phi, sides = 1, method = "recursive")
        if (length(theta) > 0)  
            x = filter(x, filter = theta, sides = 1, method = "convolution")
        x
    }
    makePulse = function(point, xtsp) {
        if (!is.ts(x)) {
            pulse = 0 * x
            pulse[point] = 1
            return(pulse)
        }
        pulse = as.numeric(0 * x)
        if (length(point) == 1) {
            x[point] = 1
            return(x)
        }
        realp = as.integer(round((point[1] + (point[2] - 1)/xtsp[3] - 
            xtsp[1]) * xtsp[3] + 1))
        pulse[realp] = 1
        pulse
    }
    upx = function(x, par) {
        ncoef = narma + ncxreg
        ind = 0
        for (modorder in transfer) {
            ind = ind + 1
            p = modorder[1]
            q = modorder[2]
            if (p > 0) 
                phi = par[ncoef + (1:p)]
            else phi = NULL
            if (q > 0) 
                theta = par[ncoef + p + (1:q)]
            else theta = NULL
            ncoef = ncoef + p + q
            x = x - armafilter(xtransf[, ind], phi = phi, theta = theta, 
                Delta = NULL)
        }
        x
    }


     btrarma = function(par, arma,trans){
      if(trans) {
      if(arma[1]) par[1:arma[1]]=levinson(par[1:arma[1]])
      if(arma[3]) par[sum(arma[c(1,2)])+seq(arma[3])]=levinson(par[sum(arma[c(1,2)])+seq(arma[3])])
}
      if(arma[1]) phi=c(1, -par[1:arma[1]]) else phi=1
      if(arma[3]) {Phi=c(1, -par[sum(arma[c(1,2)])+seq(arma[3])])
             zero=matrix(0, nrow=arma[5]-1, ncol=length(Phi))
             Phi=rbind(Phi, zero)
             Phi=as.vector(Phi)
             Phi=Phi[-seq(length(Phi)-arma[5]+2, length(Phi))]} else Phi=1
      phi=conv(phi, Phi)
      phi=-phi[-1]

      if(arma[2]) theta=c(1, par[arma[1]+seq(arma[2])]) else theta=1
      if(arma[4]) {Theta=c(1, par[sum(arma[c(1,2,3)])+seq(arma[4])])
             zero=matrix(0, nrow=arma[5]-1, ncol=length(Theta))
             Theta=rbind(Theta, zero)
             Theta=as.vector(Theta)
             Theta=Theta[-seq(length(Theta)-arma[5]+2, length(Theta))]} else Theta=1
      theta=conv(theta, Theta)
      theta=theta[-1]

      list(phi=phi, theta=theta, par=par)
      }

levinson=function(u){
#
# compute u to phi
#
#pacf=(1-exp(-u))/(1+exp(-u))
pacf=tanh(u)
p=length(u)
if (p<=1) return(pacf) else {
oldphi=pacf
for (k in 2:p) {
mult=pacf[k]
newphi=rep(NA,k)
newphi[k]=pacf[k]
for (j in 1: (k-1)) {
newphi[j]=oldphi[j]-mult*oldphi[k-j]
}
oldphi=newphi
}
}
newphi
}




     armafn <- function(p, trans) {
        par <- coef
        par[mask] <- p
        trarma <- btrarma(par, arma,trans)
        par=trarma$par

        Z <- makeARIMA(trarma[[1]], trarma[[2]], Delta, kappa)
#        r=max(length(Z$phi), length(Z$theta)+1)
#        if (r > 1) {
#            T1=Z$T[1:r,1:r]
#            V1=Z$V[1:r,1:r]
#            Z$Pn[1:r, 1:r] <- solve(diag(r*r)-T1%x%T1, as.vector(V1))
#            }
#        else if (length(Z$phi) > 0) 
#            Z$Pn[1, 1] <- 1/(1 - Z$phi^2)
#        else Z$Pn[1, 1] <- 1
#        Z$a[] <- 0
         if (ncxreg > 0) {
            nxreg = checkIO(xreg, phi = Z$phi, theta = Z$theta, 
                Delta = Z$Delta)
            x <- x - nxreg %*% par[narma + (1:ncxreg)]
        }
        if (ntransf.term) 
            x = upx(x, par)
        res=stats:::arima(x, order=order, seasonal=seasonal, fixed=par[1:narma],method='ML', include.mean=FALSE,transform.pars =FALSE)
        (-2*res$loglik-n.used*(1+log(2*pi)))/(2*n.used)
    }

     armafn1 <- function(p, trans) {
        par <- coef
        par[mask] <- p
        trarma <- btrarma(par, arma, trans)
        par=trarma$par
        Z <- makeARIMA(trarma[[1]], trarma[[2]], Delta, kappa)
#        r=max(length(Z$phi), length(Z$theta)+1)
#        if (r > 1) {
#            T1=Z$T[1:r,1:r]
#            V1=Z$V[1:r,1:r]
#            Z$Pn[1:r, 1:r] <- solve(diag(r*r)-T1%x%T1, as.vector(V1))
#            }
#        else if (length(Z$phi) > 0) 
#            Z$Pn[1, 1] <- 1/(1 - Z$phi^2)
#        else Z$Pn[1, 1] <- 1
#        Z$a[] <- 0
         if (ncxreg > 0) {
            nxreg = checkIO(xreg, phi = Z$phi, theta = Z$theta, 
                Delta = Z$Delta)
            x <- x - nxreg %*% par[narma + (1:ncxreg)]
        }
        if (ntransf.term) 
            x = upx(x, par)
        res=stats:::arima(x, order=order, seasonal=seasonal, fixed=par[1:narma],method='ML', include.mean=FALSE,transform.pars =FALSE)
        res
    }


    armaCSS <- function(p, trans) {
        par <- as.double(fixed)
        par[mask] <- p
        trarma <- btrarma(par, arma,trans)
        par=trarma$par
        Z <- makeARIMA(trarma[[1]], trarma[[2]], Delta, kappa)
                if (ncxreg > 0) {
            nxreg = checkIO(xreg, phi = trarma[[1]], theta = trarma[[2]], 
                Delta = Delta)
            x <- x - nxreg %*% par[narma + (1:ncxreg)]
        }
        if (ntransf.term) 
            x = upx(x, par)
        res=stats:::arima(x, order=order, seasonal=seasonal, fixed=par[1:narma],method='CSS', include.mean=FALSE,transform.pars =FALSE)$sigma2
        0.5*log(res)   
 }
    armaCSS1 <- function(p,trans) {
        par <- as.double(fixed)
        par[mask] <- p
        trarma <- btrarma(par, arma, trans)
        par=trarma$par
        Z <- makeARIMA(trarma[[1]], trarma[[2]], Delta, kappa)
                if (ncxreg > 0) {
            nxreg = checkIO(xreg, phi = trarma[[1]], theta = trarma[[2]], 
                Delta = Delta)
            x <- x - nxreg %*% par[narma + (1:ncxreg)]
        }
        if (ntransf.term) 
            x = upx(x, par)
        stats:::arima(x, order=order, seasonal=seasonal, fixed=par[1:narma],method='CSS', include.mean=FALSE,transform.pars =FALSE)
            }


    arCheck <- function(ar) {
        p <- max(which(c(1, -ar) != 0)) - 1
        if (!p) 
            return(TRUE)
        all(Mod(polyroot(c(1, -ar[1:p]))) > 1)
    }
    maInvert <- function(ma) {
        q <- length(ma)
        q0 <- max(which(c(1, ma) != 0)) - 1
        if (!q0) 
            return(ma)
        roots <- polyroot(c(1, ma[1:q0]))
        ind <- Mod(roots) < 1
        if (all(!ind)) 
            return(ma)
        if (q0 == 1) 
            return(c(1/ma[1], rep(0, q - q0)))
        roots[ind] <- 1/roots[ind]
        x <- 1
        for (r in roots) x <- c(x, 0) - c(0, x)/r
        c(Re(x[-1]), rep(0, q - q0))
    }

   if(missing(io) & missing(xtransf) & missing(transfer)) {
      res=stats:::arima(x=x, order=order, seasonal=seasonal, xreg=xreg, 
include.mean=include.mean, transform.pars=transform.pars,
fixed=fixed, init=init, method=method, n.cond=if(!missing(n.cond)) n.cond, optim.control=optim.control,
kappa=kappa)
      res$call=match.call()
      res$aic=res$aic-2
      return(res)
}
    
    series <- deparse(substitute(x))
    if (NCOL(x) > 1) 
        stop("only implemented for univariate time series")
    method <- match.arg(method)
    x <- as.ts(x)
    if (!is.numeric(x)) 
        stop("'x' must be numeric")
    storage.mode(x) <- "double"
    dim(x) <- NULL
    n <- length(x)
    if (!missing(order)) 
        if (!is.numeric(order) || length(order) != 3 || any(order < 
            0)) 
            stop("'order' must be a non-negative numeric vector of length 3")
    if (!missing(seasonal)) 
        if (is.list(seasonal)) {
            if (is.null(seasonal$order)) 
                stop("'seasonal' must be a list with component 'order'")
            if (!is.numeric(seasonal$order) || length(seasonal$order) != 
                3 || any(seasonal$order < 0)) 
                stop("'seasonal$order' must be a non-negative numeric vector of length 3")
        }
        else if (is.numeric(order)) {
            if (length(order) == 3) 
                seasonal <- list(order = seasonal)
            else ("'seasonal' is of the wrong length")
        }
        else stop("'seasonal' must be a list with component 'order'")
    if (is.null(seasonal$period) || is.na(seasonal$period) || 
        seasonal$period == 0) 
        seasonal$period <- frequency(x)

    arma <- as.integer(c(order[-2], seasonal$order[-2], seasonal$period, 
        order[2], seasonal$order[2]))
    narma <- sum(arma[1:4])
    xtsp <- tsp(x)
    tsp(x) <- NULL

    Pm = NULL
    if (!missing(io)) {
        Pm = sapply(io, makePulse, xtsp)
        name1 = substitute(io)
        namev = NULL
        for (i in 2:length(name1)) namev = c(namev, paste("IO", 
            deparse(name1[[i]]), sep = "-"))
        colnames(Pm) = namev
    }

    Delta <- 1
    for (i in seq(len = order[2])) Delta <- conv1(Delta)
    for (i in seq(len = seasonal$order[2])) Delta <- conv1(Delta,period=seasonal$period)
    Delta <- -Delta[-1]
    nd <- order[2] + seasonal$order[2]
    n.used <- sum(!is.na(x)) - length(Delta)

    if(!is.null(xreg)) xreg=data.frame(xreg)
    if (!is.null(xreg)) 
        xreg = as.matrix(xreg)
    xreg = cbind(xreg, Pm)

    if (is.null(xreg)) {
        ncxreg <- 0
    }
    else {
        nmxreg <- deparse(substitute(xreg))
        if (NROW(xreg) != n) 
            stop("lengths of 'x' and 'xreg' do not match")
        ncxreg <- NCOL(xreg)
        xreg <- as.matrix(xreg)
        storage.mode(xreg) <- "double"
    }
    class(xreg) <- NULL
    if (ncxreg > 0 && is.null(colnames(xreg))) 
        colnames(xreg) <- if (ncxreg == 1) 
            nmxreg
        else paste(nmxreg, 1:ncxreg, sep = "")
    if (include.mean && (nd == 0)) {
        xreg <- cbind(intercept = rep(1, n), xreg = xreg)
        ncxreg <- ncxreg + 1
    }

    ntransf = 0
    ntransf.term = 0
    ct = NULL
    if (!is.null(transfer)) {
        xtransf = as.matrix(xtransf)
        storage.mode(xtransf) = "double"
        ntransf.term = length(transfer)
        ntransf = sum(sapply(transfer, sum))
        if (missing(xtransf) && length(transfer) != ncol(xtransf)) 
            stop("no of columns of xtransf need to equal the number of transfer sub-models")
        ncoef = narma + ncxreg
        ind = 0
        if (!is.null(colnames(xtransf))) 
            tnames = colnames(xtransf)
        else tnames = paste("T", 1:ntransf.term, sep = "")
        for (modorder in transfer) {
            ind = ind + 1
            p = modorder[1]
            q = modorder[2]
            if (p > 0) 
                ct = c(ct, paste(paste(tnames[ind], "AR", sep = "-"), 
                  1:p, sep = ""))
            if (q > 0) 
                ct = c(ct, paste(paste(tnames[ind], "MA", sep = "-"), 
                  1:q - 1, sep = ""))
        }
    }
    if (method == "CSS-ML") {
        anyna <- any(is.na(x))
        if (ncxreg) 
            anyna <- anyna || any(is.na(xreg))
        if (anyna) 
            method <- "ML"
    }

    if (method == "CSS" || method == "CSS-ML") {
        ncond <- order[2] + seasonal$order[2] * seasonal$period
        ncond1 <- order[1] + seasonal$period * seasonal$order[1]
        ncond <- if (!missing(n.cond)) 
            ncond + max(n.cond, ncond1)
        else ncond + ncond1
    }
    else ncond <- 0

    if (is.null(fixed)) 
        fixed <- rep(as.numeric(NA), narma + ncxreg + ntransf)
    else if (length(fixed) != narma + ncxreg + ntransf) 
        stop("wrong length for 'fixed'")
 
    mask <- is.na(fixed)
    no.optim <- !any(mask)
    if (no.optim) 
        transform.pars <- FALSE
    if (transform.pars) {
        ind <- arma[1] + arma[2] + seq(length = arma[3])
        if (any(!mask[seq(length = arma[1])]) || any(!mask[ind])) {
            warning("some AR parameters were fixed: setting transform.pars = FALSE")
            transform.pars <- FALSE
        }
    }

    init0 <- rep(0, narma)
    parscale <- rep(1, narma)
 
    if (ncxreg) {
        cn <- colnames(xreg)
        orig.xreg <- (ncxreg == 1) || any(!mask[narma + 1:ncxreg])
        if (!orig.xreg) {
            S <- svd(na.omit(xreg))
            xreg <- xreg %*% S$v
        }
        if (include.mean && (nd == 0)) {
            fit <- lm(x ~ xreg - 1, na.action = na.omit)
            n.used <- sum(!is.na(resid(fit))) - length(Delta)
            init0 <- c(init0, coef(fit))
            ses <- summary(fit)$coefficients[, 2]
        }
        else {
            fit <- lm(x ~ xreg, na.action = na.omit)
            n.used <- sum(!is.na(resid(fit))) - length(Delta)
            init0 <- c(init0, coef(fit)[-1])
            ses <- summary(fit)$coef[-1, 2]
        }
        parscale <- c(parscale, 10 * ses)
    }

    if (ntransf) {
        init0 = c(init0, rep(0, ntransf))
    }
    if (ntransf.term) {
        for (i in 1:ntransf.term) {
            fit = lm(x ~ xtransf[, i], na.action = na.omit)
            ses = summary(fit)$coef[2, 2]
            parscale = c(parscale, rep(10 * ses, sum(transfer[[i]])))
        }
    }

    if (n.used <= 0) 
        stop("too few non-missing observations")

    if (!is.null(init)) {
        if (length(init) != length(init0)) 
            stop("'init' is of the wrong length")
        if (any(ind <- is.na(init))) 
            init[ind] <- init0[ind]
        if (method == "ML") {
            if (arma[1] > 0) 
                if (!arCheck(init[1:arma[1]])) 
                  stop("non-stationary AR part")
            if (arma[3] > 0) 
                if (!arCheck(init[sum(arma[1:2]) + 1:arma[3]])) 
                  stop("non-stationary seasonal AR part")
            if (transform.pars) init=transformpar(init, arma)
        }
    }
    else init <- init0

    coef <- as.double(fixed)

    if (!("parscale" %in% names(optim.control))) 
        optim.control$parscale <- parscale[mask]

    if (method == "CSS") {
        res <- if (no.optim) 
            list(convergence = 0, par = numeric(0), value = armaCSS(numeric(0)))
        else optim(init[mask], armaCSS, method = "BFGS", hessian = TRUE, 
            control = optim.control, trans=transform.pars)
        if (res$convergence > 0) 
            warning("possible convergence problem: optim gave code=", 
                res$convergence)
        coef[mask] <- res$par
        val <- armaCSS1(res$par, transform.pars)
        sigma2 <- val$sigma2
        mod = val$mod  
        var <- if (no.optim) 
            numeric(0)
        else solve(res$hessian * n.used)
    }
    else {
        if (method == "CSS-ML") {
            res <- if (no.optim) 
                list(convergence = 0, par = numeric(0), value = armaCSS(numeric(0), trans=FALSE))
            else optim(init[mask], armaCSS, method = "BFGS", 
                hessian = FALSE, control = optim.control, trans=transform.pars)
           if (res$convergence == 0) { 
                if(transform.pars) res$par=undotrans(res$par, arma)
                init[mask] <- res$par} 
            if (arma[1] > 0) 
                if (!arCheck(init[1:arma[1]])) 
                  stop("non-stationary AR part from CSS")
            if (arma[3] > 0) 
                if (!arCheck(init[sum(arma[1:2]) + 1:arma[3]])) 
                  stop("non-stationary seasonal AR part from CSS")
            ncond <- 0
        }

         if (transform.pars) {
             init=transformpar(init, arma)
            if (arma[2] > 0) {
                ind <- arma[1] + 1:arma[2]
                init[ind] <- maInvert(init[ind])
            }
            if (arma[4] > 0) {
                ind <- sum(arma[1:3]) + 1:arma[4]
                init[ind] <- maInvert(init[ind])
            }
        }
        trarma <- btrarma(init, arma, transform.pars)
        mod <- makeARIMA(trarma[[1]], trarma[[2]], Delta, kappa)
        res <- if (no.optim) 
            list(convergence = 0, par = numeric(0), value = armafn(numeric(0), 
                as.logical(transform.pars)))
        else optim(init[mask], armafn, method = "BFGS", hessian = TRUE, 
            control = optim.control, trans = as.logical(transform.pars))
        
        if (res$convergence > 0) 
            warning("possible convergence problem: optim gave code=", 
                res$convergence)
       coef[mask] <- res$par

        if (transform.pars) {
            if (arma[2] > 0) {
                ind <- arma[1] + 1:arma[2]
                if (all(mask[ind])) 
                  coef[ind] <- maInvert(coef[ind])
            }
            if (arma[4] > 0) {
                ind <- sum(arma[1:3]) + 1:arma[4]
                if (all(mask[ind])) 
                  coef[ind] <- maInvert(coef[ind])
            }
            if (any(coef[mask] != res$par)) {
                oldcode <- res$convergence
                res <- optim(coef[mask], armafn, method = "BFGS", 
                  hessian = TRUE, control = list(maxit = 0, parscale = optim.control$parscale), 
                  trans = TRUE)
                res$convergence <- oldcode
                coef[mask] <- res$par
            }

             A = difpar(as.double(coef), arma)
             A <- A[mask, mask]
            var <- t(A) %*% solve(res$hessian * n.used) %*% A
            coef = undotrans(coef, arma)
        }
        else var <- if (no.optim) 
            numeric(0)
        else solve(res$hessian * n.used)
         val = armafn1(res$par,as.logical(transform.pars)) 
         sigma2 <- val$sigma2 
         mod = val$mod   
}

    value <-  2 * n.used * res$value + n.used + n.used * log(2 * 
        pi)    
    aic <- if (method != "CSS") 
        value + 2 * sum(mask)
    else NA
    nm <- NULL
    if (arma[1] > 0) 
        nm <- c(nm, paste("ar", 1:arma[1], sep = ""))
    if (arma[2] > 0) 
        nm <- c(nm, paste("ma", 1:arma[2], sep = ""))
    if (arma[3] > 0) 
        nm <- c(nm, paste("sar", 1:arma[3], sep = ""))
    if (arma[4] > 0) 
        nm <- c(nm, paste("sma", 1:arma[4], sep = ""))
    if (ncxreg > 0) {
        nm <- c(nm, cn)
        if (!orig.xreg) {
            ind <- narma + 1:ncxreg
            coef[ind] <- S$v %*% coef[ind]
            A <- diag(length(coef))
            A[ind, ind] <- S$v
            A <- A[mask, mask]
            var <- A %*% var %*% t(A)
        }
    }
    nm = c(nm, ct)
    names(coef) <- nm


    if (!no.optim) 
        dimnames(var) <- list(nm[mask], nm[mask])
    resid <- val$residuals
    tsp(resid) <- xtsp
    class(resid) <- "ts"
    res <- list(coef = coef, sigma2 = sigma2, var.coef = var, 
        mask = mask, loglik = -0.5 * value, aic = aic, arma = arma, 
        residuals = resid, call = match.call(), series = series, 
        code = res$convergence, n.cond = ncond, model = mod)
    class(res) <- c("Arimax", "Arima")
    res
}
