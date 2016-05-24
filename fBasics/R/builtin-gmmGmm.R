
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################
# FUNCTION:             
#  .gmm
#  .HAC
#  .weightsAndrews2
#  .bwAndrews2
#  .summary.gmm
#  .lintest
################################################################################


# Code borrowed from 
#   R's contributed package "gmm" written by Pierre Chausse.


# Rmetrics:
#   Note that gmm is not available on Debian as of 2009-04-28. 
#   To run these functions under Debian/Rmetrics we have them    
#   implemented here as a builtin.
#   We also made modifications for tailored usage with Rmetrics. 
#   Note that the dependences in the original package requires zoo
#   which may create conflicts with Rmetrics timeDate/timeSeries.


# Package: gmm
# Version: 1.0-4
# Date: 2009-03-10
# Title: Generalized Method of Moments and Generalized Empirical Likelihood
# Author: Pierre Chausse <pierre.chausse@uqam.ca>
# Maintainer: Pierre Chausse <pierre.chausse@uqam.ca>
# Description: It is a complete suite to estimate models based on moment 
#     conditions. It includes the  two step Generalized method of moments 
#     (GMM) of Hansen(1982), the iterated GMM and continuous updated 
#     estimator (CUE) of Hansen-Eaton-Yaron(1996) and several methods that 
#     belong to the Generalized Empirical Likelihood (GEL) family of estimators, 
#     as presented by Smith(1997), Kitamura(1997), Newey-Smith(2004) and 
#     Anatolyev(2005). 
# Depends: R (>= 2.0.0), sandwich, tseries, mvtnorm
# Imports: stats
# License: GPL (>= 2)


# ------------------------------------------------------------------------------


.gmm <- 
function(g, x, t0=NULL, gradv=NULL, type=c("twoStep","cue","iterative"), 
    wmatrix = c("optimal","ident"),  vcov=c("HAC","iid"), 
    kernel=c("Quadratic Spectral","Truncated", "Bartlett", "Parzen", 
    "Tukey-Hanning"),crit=10e-7,bw = .bwAndrews2, 
    prewhite = FALSE, ar.method = "ols", approx="AR(1)",tol = 1e-7, 
    itermax=100,intercept=TRUE,optfct=c("optim","optimize"), ...)
{
    type <- match.arg(type)
    kernel <- match.arg(kernel)
    vcov <- match.arg(vcov)
    wmatrix <- match.arg(wmatrix)
    optfct <- match.arg(optfct)
    typeg=0

    if (is(g,"formula"))
        {
        typeg=1
        dat <- .get_dat(g,x,intercept=intercept)
        x <- dat$x

        g <- function(tet,x,ny=dat$ny,nh=dat$nh,k=dat$k)
            {
            tet <- matrix(tet,ncol=k)
            if (intercept)
                {
                e <- x[,1:ny] -  x[,(ny+1):(ny+k)]%*%t(tet)
                gt <- e
                for (i in 2:nh)
                    {
                    gt <- cbind(gt,e*x[,(ny+k+i)])
                    }
                }
            if (!intercept)
                e <- x[,1:ny] -  x[,(ny+1):(ny+k)]%*%t(tet)
                gt <- e*x[,ny+k+1]
                if (nh > 1)
                    {    
                    for (i in 2:nh)
                        {
                        gt <- cbind(gt,e*x[,(ny+k+i)])
                        }
                    }
            return(gt)
            }
        gradv <- function(x,ny=dat$ny,nh=dat$nh,k=dat$k)
            {
            dgb <- -(t(x[,(ny+k+1):(ny+k+nh)])%*%x[,(ny+1):(ny+k)])%x%
                diag(rep(1,ny))/nrow(x)
            return(dgb)
            }
        tetlin <- function(x,w,ny=dat$ny,nh=dat$nh,k=dat$k)
            {
            n <- nrow(x)
            ym <- as.matrix(x[,1:ny])
            xm <- as.matrix(x[,(ny+1):(ny+k)])
            hm <- as.matrix(x[,(ny+k+1):(ny+k+nh)])
            whx <- solve(w,(crossprod(hm,xm)%x%diag(ny)))    
            wvecyh <- solve(w,matrix(crossprod(ym,hm),ncol=1))
            dg <- gradv(x)
            xx <- crossprod(dg,whx)
            par <- solve(xx,crossprod(dg,wvecyh))
            gb <- matrix(colSums(g(par,x))/n,ncol=1)
            value <- crossprod(gb,solve(w,gb)) 
            res <- list(par=par,value=value)
            return(res)
            }
        }
if (optfct == "optimize")
    {
    n = nrow(g(t0[1],x))
    q = ncol(g(t0[1],x))
    k = 1
    k2 <- k
    df <- q-k
    }
else
    {
    if (typeg)
        {
        k <- dat$k
        k2 <- k*dat$ny
        n <- nrow(x)
        q <- dat$ny*dat$nh
        df <- q-k*dat$ny
        }
    else
        {
        n = nrow(g(t0,x))
        q = ncol(g(t0,x))
        k = length(t0)
        k2 <- k
        df <- q-k
        }
    }
obj1 <- function(thet)
    {
    gt <- g(thet,x)
    gbar <- as.vector(colMeans(gt))
    obj <- crossprod(gbar,solve(w,gbar))
    return(obj)
    }
iid <- function(thet)
    {
    gt <- g(thet,x)
    v <- crossprod(gt,gt)/n
    return(v)
    }
Gf <- function(thet)
    {
    myenv <- new.env()
    assign('x',x,envir=myenv)
    assign('thet',thet,envir=myenv)
    barg <- function(x,thet)
        {
        gt <- g(thet,x)
        gbar <- as.vector(colMeans(gt))
        return(gbar)
        }
    G <- attr(numericDeriv(quote(barg(x,thet)),"thet",myenv),"gradient")
    return(G)
    }

if (q == k2 | wmatrix == "ident")
    {
    if (typeg)
        {
        w <- diag(q)
        res <- tetlin(x,w)
        z = list(par=res$par,objective=res$value)
        }
    else
        {
        w=diag(rep(1,q))
        if (optfct == "optim")
            res <- optim(t0,obj1, ...)
        else
            {
            res <- optimize(obj1,t0, ...)
            res$par <- res$minimum
            res$value <- res$objective
            }    
        z = list(par=res$par,objective=res$value)    
        }
    }
else
    {
    if (type=="twoStep")
        {
        w=diag(rep(1,q))
        if (typeg)
            {
            res1 <- tetlin(x,w)
            }
        else
            {
            if (optfct == "optim")
                res1 <- optim(t0,obj1, ...)
            else
                {
                res1 <- optimize(obj1,t0, ...)
                res1$par <- res1$minimum
                res1$value <- res1$objective
                }    
            }
        if (vcov == "iid")
            w <- iid(res1$par)
        if (vcov == "HAC")
            w <- .HAC(g(res1$par,x), kernel=kernel, bw=bw,
                prewhite=prewhite,ar.method=ar.method,approx=approx,tol=tol)
        if (typeg)
            {
            res2 <- tetlin(x,w)
            }
        else
            {
            if (optfct == "optim")
                res2 <- optim(res1$par,obj1, ...)
            else
                {
                res2 <- optimize(obj1,t0, ...)
                res2$par <- res2$minimum
                res2$value <- res2$objective
                }    
            }
        z = list(par=res2$par,objective=res2$value)    
        }
    if (type=="cue")
        {
        obj_cue <- function(thet)
            {
            gt <- g(thet,x)
            gbar <- as.vector(colMeans(gt))
            if (vcov == "iid")
                w2 <- iid(thet)
            if (vcov == "HAC")
                w2 <- .HAC(g(thet,x), kernel=kernel, bw=bw,
                    prewhite=prewhite,ar.method=ar.method,approx=approx,tol=tol)
            obj <- crossprod(gbar,solve(w2,gbar))
            return(obj)
            }    
        if (typeg)
            {
            if (is.null(t0))
                t0 <- tetlin(x,diag(rep(1,q)))$par
            if (optfct == "optim")
                {
                res2 <- optim(t0,obj_cue, ...)
                }    
            else
                {
                res2 <- optimize(obj_cue,t0, ...)
                res2$par <- res2$minimum
                res2$value <- res2$objective
                }
            }
        else
            {
            if (optfct == "optim")
                res2 <- optim(t0,obj_cue, ...)
            else
                {
                res2 <- optimize(obj_cue,t0, ...)
                res2$par <- res2$minimum
                res2$value <- res2$objective
                }    
            }
        z = list(par=res2$par,objective=res2$value)    
        }
    if (type=="iterative")
        {
        w=diag(rep(1,q))
        if (typeg)
            res <- tetlin(x,w)
        else
            {
            if (optfct == "optim")
                res <- optim(t0,obj1, ...)
            else
                {
                res <- optimize(obj1,t0, ...)
                res$par <- res$minimum
                res$value <- res$objective
                }    
            }
        ch <- 100000
        j <- 1
        while(ch>crit)
        {
            tet <- res$par
            if (vcov == "iid")
                w2 <- iid(tet)
            if (vcov == "HAC")
                w2 <- .HAC(g(tet,x), kernel=kernel, bw=bw,
                    prewhite=prewhite,ar.method=ar.method,approx=approx,tol=tol)
            if (typeg)
                res <- tetlin(x,w2)
            else
                {    
                if (optfct == "optim")
                    res <- optim(tet,obj1, ...)
                else
                    {
                    res <- optimize(obj1,t0, ...)
                    res$par <- res$minimum
                    res$value <- res$objective
                    }    
                }
            ch <- crossprod(abs(tet-res$par)/tet,abs(tet-res$par)/tet)
            if (j>itermax)
                {
                cat("No convergence after ",itermax," iterations")
                ch <- crit
                }
            j <- j+1    
        }

        z = list(par=res$par,objective=res$value)    
        }
    }

    if (!is.function(gradv)) 
        G <- Gf(z$par)
    else
        if (typeg)
            G <- gradv(x)
        else    
            G <- gradv(z$par,x)

    if (vcov == "iid")
        v <- iid(z$par)/n
    else
        v <- .HAC(g(z$par,x), kernel=kernel, bw=bw,prewhite=prewhite,
            ar.method=ar.method,approx=approx,tol=tol)/n
    
    if (wmatrix == "optimal")
        {
        z$vcov <- solve(crossprod(G,solve(v,G)))
        }
    else
        {
        GGG <- solve(crossprod(G),t(G))
        z$vcov <- GGG%*%v%*%t(GGG)
        }

z$gt <- g(z$par,x)
if (typeg==0)
    names(z$par) <- paste("Theta[",1:k,"]",sep="")
if (typeg == 1)
    {
    namex <- colnames(dat$x[,(dat$ny+1):(dat$ny+dat$k)])
    nameh <- colnames(dat$x[,(dat$ny+dat$k+1):(dat$ny+dat$k+dat$nh)])
    if (dat$ny > 1)
        {
        namey <- colnames(dat$x[,1:dat$ny])
        names(z$par) <- paste(rep(namey,dat$k),"_",rep(namex,rep(dat$ny,dat$k)),sep="")
        colnames(z$gt) <- paste(rep(namey,dat$nh),"_",rep(nameh,rep(dat$ny,dat$nh)),sep="")
        }
    if (dat$ny == 1)
        {
        names(z$par) <- namex
        colnames(z$gt) <- nameh
        }
    }

    dimnames(z$vcov) <- list(names(z$par),names(z$par))
    z$df <- df
    z$k <- k
    z$n <- n
    z$met <- type
    z$kernel <- kernel
    class(z) <- "gmm"
    z
}


# ------------------------------------------------------------------------------


.HAC <- 
function(x, weights = .weightsAndrews2, bw = .bwAndrews2, 
    prewhite = FALSE, ar.method = "ols", kernel=c("Quadratic Spectral", 
    "Truncated", "Bartlett", "Parzen", "Tukey-Hanning"), approx="AR(1)",
    tol = 1e-7) 
{
    n.orig <- n <- nrow(x)
    k <- ncol(x)
    kernel=match.arg(kernel)    
    if(prewhite > 0) 
    {
        var.fit <- ar(x, order.max = prewhite, demean = FALSE, aic = FALSE, 
            method = ar.method)
        if(k > 1) D <- solve(diag(ncol(x)) - apply(var.fit$ar, 2:3, sum))
         else D <- as.matrix(1/(1 - sum(var.fit$ar)))
    x <- as.matrix(na.omit(var.fit$resid))
    n <- n - prewhite
    }
    weights <- weights(x, ar.method = ar.method,kernel=kernel,bw=bw, 
        approx = approx, prewhite = 1, tol = tol)
    if (length(weights) > n) 
    {
        warning("more weights than observations, only first n used")
        weights <- weights[1:n]
    }
    utu <- 0.5 * crossprod(x) * weights[1]
    wsum <- n * weights[1]/2
    w2sum <- n * weights[1]^2/2
    if (length(weights) > 1) {
        for (ii in 2:length(weights)) {
            utu <- utu + weights[ii] * crossprod(x[1:(n - 
                ii + 1), , drop = FALSE], x[ii:n, , drop = FALSE])
            wsum <- wsum + (n - ii + 1) * weights[ii]
            w2sum <- w2sum + (n - ii + 1) * weights[ii]^2
        }
    }
    utu <- utu + t(utu)
    
    if(prewhite > 0) {
    utu <- crossprod(t(D), utu) %*% t(D)
     }
    wsum <- 2 * wsum
    w2sum <- 2 * w2sum
    bc <- n^2/(n^2 - wsum)
    df <- n^2/w2sum
    rval <- utu/n.orig
    
    return(rval)
}


# ------------------------------------------------------------------------------


.weightsAndrews2 <- 
function(x, bw = .bwAndrews2, kernel = c("Quadratic Spectral", 
    "Truncated", "Bartlett", "Parzen", "Tukey-Hanning"), approx = c("AR(1)", 
    "ARMA(1,1)"), prewhite = 1, ar.method = "ols", tol = 1e-7, verbose = FALSE)
{
    kernel <- match.arg(kernel)
    approx=match.arg(approx)

    if (is.function(bw)) 
        bw <- bw(x, kernel = kernel, prewhite = prewhite, 
            ar.method = ar.method, approx=approx)
    n <- NROW(x) 
    weights <- .kweights(0:(n - 1)/bw, kernel = kernel)
    weights <- weights[1:max(which(abs(weights) > tol))]
    return(weights)
}


# ------------------------------------------------------------------------------


.bwAndrews2 <- 
function(x, kernel = c("Quadratic Spectral", 
    "Truncated", "Bartlett", "Parzen", "Tukey-Hanning"), approx = c("AR(1)", 
    "ARMA(1,1)"), prewhite = 1, ar.method = "ols") 
{
    kernel <- match.arg(kernel)
    approx <- match.arg(approx)
    n <- nrow(x)
    k <- ncol(x)

    if (approx == "AR(1)") {
        fitAR1 <- function(x) {
            rval <- ar(x, order.max = 1, aic = FALSE, method = "ols")
            rval <- c(rval$ar, sqrt(rval$var.pred))
            names(rval) <- c("rho", "sigma")
            return(rval)
        }
        ar.coef <- apply(x, 2, fitAR1)
        denum <- sum((ar.coef["sigma", ]/(1 - ar.coef["rho", 
            ]))^4)
        alpha2 <- sum(4 * ar.coef["rho", ]^2 * ar.coef["sigma", 
            ]^4/(1 - ar.coef["rho", ])^8)/denum
        alpha1 <- sum(4 * ar.coef["rho", ]^2 * ar.coef["sigma", 
            ]^4/((1 - ar.coef["rho", ])^6 * (1 + ar.coef["rho", 
            ])^2))/denum
    }
    else {
        fitARMA11 <- function(x) {
            rval <- arima(x, order = c(1, 0, 1), include.mean = FALSE)
            rval <- c(rval$coef, sqrt(rval$sigma2))
            names(rval) <- c("rho", "psi", "sigma")
            return(rval)
        }
        arma.coef <- apply(x, 2, fitARMA11)
        denum <- sum(((1 + arma.coef["psi", ]) * arma.coef["sigma", 
            ]/(1 - arma.coef["rho", ]))^4)
        alpha2 <- sum(4 * ((1 + arma.coef["rho", ] * 
            arma.coef["psi", ]) * (arma.coef["rho", ] + arma.coef["psi", 
            ]))^2 * arma.coef["sigma", ]^4/(1 - arma.coef["rho", 
            ])^8)/denum
        alpha1 <- sum(4 * ((1 + arma.coef["rho", ] * 
            arma.coef["psi", ]) * (arma.coef["rho", ] + arma.coef["psi", 
            ]))^2 * arma.coef["sigma", ]^4/((1 - arma.coef["rho", 
            ])^6 * (1 + arma.coef["rho", ])^2))/denum
    }
    rval <- switch(kernel, Truncated = {
        0.6611 * (n * alpha2)^(1/5)
    }, Bartlett = {
        1.1447 * (n * alpha1)^(1/3)
    }, Parzen = {
        2.6614 * (n * alpha2)^(1/5)
    }, "Tukey-Hanning" = {
        1.7462 * (n * alpha2)^(1/5)
    }, "Quadratic Spectral" = {
        1.3221 * (n * alpha2)^(1/5)
    })
   return(rval)
}


# ------------------------------------------------------------------------------


.summary.gmm <- 
function(object, interval=FALSE, ...)
    {
    z <- object
    se <- sqrt(diag(z$vcov))
    par <- z$par
    tval <- par/se
    j <- z$objective*z$n
    ans <- list(met=z$met,kernel=z$kernel,algo=z$algo)
    names(ans$met) <- "GMM method"
    names(ans$kernel) <- "kernel for cov matrix"
        
    ans$par <- round(cbind(par,se, tval, 2 * pnorm(abs(tval), lower.tail = FALSE)),5)

        dimnames(ans$par) <- list(names(z$par), 
        c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))

    ans$J_test <- noquote(paste("Test-J degrees of freedom is ",z$df,sep=""))
    ans$j <- noquote(cbind(j,ifelse(z$df>0,pchisq(j,z$df,lower.tail = FALSE),"*******")))
    dimnames(ans$j) <- list("Test E(g)=0:  ",c("J-test","Pz(>j)"))
    
    if (interval != FALSE)
        {
        zs <- qnorm((1-interval)/2,lower.tail=FALSE)
        ch <- zs*se
        ans$interval <- cbind(par-ch,par+ch)
        dimnames(ans$interval) <- list(names(par),c("Theta_lower","Theta_upper"))
        }
    class(ans) <- "summary.gmm"
    
    ans
}


# ------------------------------------------------------------------------------
    
    
.lintest <- 
function(object,R,c)
{
    z <- object
    dl <- nrow(R)
    rest <- R%*%z$par - c
    if (class(z)=="gmm")
        vcov_par <- z$vcov
    else
        vcov_par <- z$vcov_par

    vcov <- R%*%vcov_par%*%t(R)
    h0 <- matrix(rep(NA,nrow(R)),ncol=1)    
    for (i in 1:nrow(R))
        {
        testnames <- names(z$par)[R[i,]!=0]
        rn <- R[i,][R[i,]!=0]
        if (rn[1] != 1)
            h0[i] <- paste(rn[1],"*",testnames[1],sep="")
        else
            h0[i] <- testnames[1]
        if (length(rn) > 1)
            {
            for (j in 2:length(rn))
                {
                if (rn[j] >= 0)
                    {    
                    if (rn[j] != 1)                    
                        h0[i] <- paste(h0[i]," + ",rn[j],"*",testnames[j],sep="")
                    else
                        h0[i] <- paste(h0[i]," + ",testnames[j],sep="")
                    }
                else
                    {    
                    if (abs(rn[j]) != 1)                    
                        h0[i] <- paste(h0[i]," - ",abs(rn[j]),"*",testnames[j],sep="")                    
                    else
                        h0[i] <- paste(h0[i]," - ",testnames[j],sep="")
                    }
                }
            }
        h0[i] <- paste(h0[i]," = ",c[i],sep="")
        }
    rh <- solve(vcov,rest)        
    ans <- list(description <- "Wald test for H0: R(Theta)=c")
    ans$H0 <- noquote(h0)
    colnames(ans$H0) <- "Null Hypothesis"        
    wt <- crossprod(rest,rh)
    pv <- pchisq(wt,dl,lower.tail=FALSE)
    res <- cbind(wt,pv)
    dimnames(res) <- list("Wald test", c("Statistics","P-Value"))
    ans$result <- res
    ans
} 


################################################################################

    