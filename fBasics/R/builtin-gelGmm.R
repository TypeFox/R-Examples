
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
#  .rho
#  .get_lamb
#  .gel
#  .smooth_g
#  .bwNeweyWest2
#  .summary.gel
#  .get_dat
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


.rho <- 
function(x, lamb, derive = 0, type = c("EL", "ET", "CUE"), drop = TRUE)
{
    type <- match.arg(type)
    lamb <- matrix(lamb,ncol=1)
    gml <- x%*%lamb
    ch <- 0
    if (derive==0)
        {
        if (type == "EL")
            {
            ch <- sum(gml>=1)
            if (drop)
                {                
                gml <- (gml<1)*gml
                rhomat <- log(1-gml) 
                }
            else
                {
                if (ch>0)
                    rhomat <- NaN
                else
                    rhomat <- log(1-gml) 
                }
            }
        if (type == "ET")
            rhomat <- -exp(gml)
        
        if (type == "CUE")
            rhomat <- -gml -0.5*gml^2
        }
    if (derive==1)
        {
        if (type == "EL")
            rhomat <- -1/(1-gml) 
            
        if (type == "ET")
            rhomat <- -exp(gml)
        
        if (type == "CUE")
            rhomat <- -1 -gml
        }
    if (derive==2)
        {
        if (type == "EL")
            rhomat <- -1/(1-gml)^2 
            
        if (type == "ET")
            rhomat <- -exp(gml)
        
        if (type == "CUE")
            rhomat <- -1
        }
    rhom <-list(ch = ch, rhomat=rhomat) 
    return(rhom)
    }


# ------------------------------------------------------------------------------

    
.get_lamb <- 
function(g, tet, x,type = c('EL', 'ET', 'CUE'), tol_lam=1e-12, 
    maxiterlam=1000, tol_obj = 1e-7)
{
    type <- match.arg(type)    
    gt <- g(tet,x)
    n <- nrow(gt)
    tol_cond=1e-12
    gb <- colMeans(gt)
    khat <- crossprod(gt)/n
    lamb0 <- -solve(khat,gb)
    conv_mes <- "Normal convergence" 

    singular <-0
    crit <-1e30
    crit0 <- crit
    dcrit <- 10
    dgblam <- -10
    gblam0 <- NULL

    j <- 1
    while ((crit > tol_lam*( 1+sqrt( crossprod(lamb0) ) ) ) & (j<=maxiterlam))
        { 
        rho2 <- as.numeric(.rho(gt,lamb0,derive=2,type=type)$rhomat)
        rho1 <- as.numeric(.rho(gt,lamb0,derive=1,type=type)$rhomat)
        gblam <- colMeans(rho1*gt)
        klam <- crossprod(rho2*gt,gt)/n
        chklam <- sum(abs(klam))
        if (!is.null(gblam0))
            dgblam <- crossprod(gblam)-crossprod(gblam0)
        
        if (is.na(chklam) | chklam == 0 | chklam == Inf |  dgblam>0 | dgblam == Inf | is.na(dgblam) | dcrit < 0)
            {
            lamb1 <- rep(sqrt(1/n),length(lamb0))
            crit <- 0
            singular=2
            conv_mes <- "The algorithm produced singular system,  NaN or Inf" 
            }
        else
            {
            if (rcond(klam)>tol_cond)
                {
                lamb1 <- lamb0-solve(klam,gblam)
                crit <- sqrt(crossprod(lamb0-lamb1))
                lamb0 <- lamb1
                }
            else
                {
                lamb1 <- rep(sqrt(1/n),length(lamb0))
                crit <- 0
                singular=2
                conv_mes <- "The algorithm produced singular system" 
                }
            }
        gblam0 <- gblam
        j <- j+1
        dcrit<- crit0-crit
        crit0 <- crit
        }
    z <- list("lambda"=lamb1,singular=singular,conv_mes=conv_mes)
    if (j>maxiterlam | max(abs(gblam))>tol_obj)
        {
        singular <- 1
        conv_mes <- "No convergence after 'maxiterlam' iterations"
        z$singular <- singular        
        }
        z$obj <- crossprod(gblam)
    return(z)
    }


# ------------------------------------------------------------------------------
    
        
.gel <- 
function(g,x,tet0,gradv=NULL,smooth=FALSE,type=c("EL","ET","CUE","ETEL"), 
    vcov=c("HAC","iid"), kernel = c("Bartlett", "Parzen", "Truncated", 
    "Tukey-Hanning"), bw = .bwAndrews2, approx = c("AR(1)", "ARMA(1,1)"), 
    prewhite = 1, ar.method = "ols", tol_weights = 1e-7, tol_lam=1e-9, 
    tol_obj = 1e-9, tol_mom = 1e-9,maxiterlam=1000, constraint = FALSE,
    intercept = TRUE, optfct = c("optim", "optimize"), 
    optlam = c("iter", "numeric"), ...)
    {
    vcov=match.arg(vcov)
    type <- match.arg(type)
    optfct <- match.arg(optfct)
    optlam <- match.arg(optlam)
    weights = .weightsAndrews2
    if (type == "ETEL")
        {
        typel <- "ET"
        typet <- "EL"    
        }
    else
        {
        typel <- type
        typet <- type
        }
    approx <- match.arg(approx)
    kernel <- match.arg(kernel)
    k <- length(tet0)

    typeg=0
    if (is(g,"formula"))
        {
        typeg=1
        dat <- .get_dat(g, x, intercept = intercept)
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
        gradv <- function(tet,x,ny=dat$ny,nh=dat$nh,k=dat$k)
            {
            tet <- matrix(tet,ncol=k)
            dgb <- -(t(x[,(ny+k+1):(ny+k+nh)])%*%x[,(ny+1):(ny+k)]) %x%
                diag(rep(1,ny))/nrow(x)
            return(dgb)
            }
        }    
        if (typeg)
            n <- nrow(x)
        else
            n = nrow(g(tet0,x))

    if (smooth)
        {
        g1 <- g
        rgmm <- .gmm(g,x,tet0,wmatrix="ident")

        if (is.function(weights))
            w <- weights(g(rgmm$par,x), kernel=kernel, bw=bw,
                prewhite = prewhite,ar.method=ar.method,approx=approx,
                tol=tol_weights)
        else
            w <- weights
        sg <- function(thet,x)
            {
            gf <- g1(thet,x)
            gt <- .smooth_g(gf, weights=w)$smoothx 
            return(gt)
            }
        g <- sg
        }    
    lll <- 1
    thetf <- function(tet)
        {
        if (optlam == "iter")
            {
            lamblist <- .get_lamb(g,tet,x,type=typel,tol_lam=tol_lam,
                maxiterlam=maxiterlam,tol_obj=tol_obj)
            lamb <- lamblist$lambda
            gt <- g(tet,x)
            pt <- -.rho(gt,lamb,type=typet,derive=1)$rhomat/nrow(gt)
            checkmom <- sum(as.numeric(pt)*gt)
            if (lamblist$singular==0)        
                p <- sum(.rho(gt,lamb,type=typet)$rhomat) + abs(checkmom)/tol_mom
            if (lamblist$singular==1)        
                p <- sum(.rho(gt,lamb,type=typet)$rhomat) + abs(checkmom)/tol_mom + lamblist$obj/tol_mom
            if (lamblist$singular==2)        
                p <- 1e50*lll
            lll <- lll+1
            }
        else
            {
            gt <- g(tet,x)
            rhofct <- function(lamb)
                {
                rhof <- -sum(.rho(gt,lamb,type=typet)$rhomat)
                return(rhof)
                }
            if (ncol(gt)>1)
                rlamb <- optim(rep(0,ncol(gt)),rhofct,...)
            else
                {
                rlamb <- optimize(rhofct,c(-1,1))
                rlamb$par <- rlamb$minimum
                rlamb$value <- rlamb$objective
                }
            lamb <- rlamb$par
            pt <- -.rho(gt,lamb,type=typet,derive=1)$rhomat/nrow(gt)
            checkmom <- sum(as.numeric(pt)*gt)
            p <- -rlamb$value + (checkmom)^2/tol_mom + (sum(as.numeric(pt))-1)^2/tol_mom
            }
        return(p)
        }

    if (!constraint)
        if (optfct == "optim")
            res <- optim(tet0,thetf,...)
        else
            {
            res <- optimize(thetf,tet0,...)
            res$par <- res$minimum
            res$convergence <- "Pas applicable"
            }
    else
        res <- constrOptim(tet0,thetf,...)


    rlamb <- .get_lamb(g,res$par,x,type=typel,tol_lam=tol_lam,maxiterlam=maxiterlam,tol_obj=tol_obj)
    z <- list(par=res$par,lambda=rlamb$lam,conv_lambda=rlamb$conv_mes,conv_par=res$convergence)
    z$foc_lambda <- rlamb$obj
    z$type <- type
    z$gt <- g(z$par,x)
    rhom <- .rho(z$gt,z$lambda,type=typet)
    z$pt <- -.rho(z$gt,z$lambda,type=typet,derive=1)$rhomat/n
    z$conv_moment <- colSums(as.numeric(z$pt)*z$gt)
    z$conv_pt <- sum(as.numeric(z$pt))
     
    z$objective <- sum(as.numeric(rhom$rhomat)-.rho(1,0,type=typet)$rhomat)/n

    if (type == "EL")    
        {
        z$badrho <- rhom$ch
        names(z$badrho) <- "Number_of_bad_rho"
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
    if (!is.function(gradv)) 
        G <- Gf(z$par)
    else
        G <- gradv(z$par,x)
    if (vcov == "iid")
        khat <- crossprod(z$gt)
    else
        khat <- .HAC(g(z$par,x), kernel=kernel, bw=bw,prewhite=prewhite,
            ar.method=ar.method,approx=approx,tol=tol_weights)

    kg <- solve(khat,G)
    z$vcov_par <- solve(crossprod(G,kg))/n
    z$vcov_lambda <- ((solve(khat)-kg%*%z$vcov_par%*%t(kg)))/n
    
    if (smooth) z$weights<-w

    if (typeg ==0)
        {
        names(z$par) <- paste("Theta[",1:k,"]",sep="")
        colnames(z$gt) <- paste("gt[",1:ncol(z$gt),"]",sep="")
        names(z$lambda) <- paste("Lambda[",1:ncol(z$gt),"]",sep="")
        }
    if (typeg == 1)
        {
        namex <- colnames(dat$x[,(dat$ny+1):(dat$ny+dat$k)])
        nameh <- colnames(dat$x[,(dat$ny+dat$k+1):(dat$ny+dat$k+dat$nh)])
        if (dat$ny > 1)
            {
            namey <- colnames(dat$x[,1:dat$ny])
            names(z$par) <- paste(rep(namey,dat$k),"_",rep(namex,rep(dat$ny,dat$k)),sep="")
            colnames(z$gt) <- paste(rep(namey,dat$nh),"_",rep(nameh,rep(dat$ny,dat$nh)),sep="")
            names(z$lambda) <- paste("Lam(",rep(namey,dat$nh),"_",rep(nameh,rep(dat$ny,dat$nh)),")",sep="")
            }
        if (dat$ny == 1)
            {
            names(z$par) <- namex
            colnames(z$gt) <- nameh
            names(z$lambda) <- nameh
            }
        }
    dimnames(z$vcov_par) <- list(names(z$par),names(z$par))
    dimnames(z$vcov_lambda) <- list(names(z$lambda),names(z$lambda))

    class(z) <- "gel"
    return(z)
    }

    
# ------------------------------------------------------------------------------


.smooth_g <- 
function(x, bw = .bwAndrews2, prewhite = 1, ar.method = "ols",
    weights=.weightsAndrews2,
    kernel=c("Bartlett","Parzen","Truncated","Tukey-Hanning"), 
    approx = c("AR(1)","ARMA(1,1)"), tol = 1e-7) 
{
    kernel <- match.arg(kernel)
    approx <- match.arg(approx)
        
    n <- nrow(x)
    if (is.function(weights))
        {
            w <- weights(x, bw = bw, kernel = kernel,  
            prewhite = prewhite, ar.method = ar.method, tol = tol, 
            verbose = FALSE, approx = approx)
        }
        else
            w <- weights


    rt <- length(w)
    if (rt >= 2)
        {
        w <- c(w[rt:2],w)
        w <- w / sum(w)
        rt <- rt-1
        sgt <- function(t) crossprod(x[(t-rt):(t+rt),],w)
        x[(rt+1):(n-rt),] <- t(sapply((rt+1):(n-rt),sgt))
        sx <- list("smoothx"=x,"kern_weights"=w)
        return(sx)        
        }
    else
        sx <- list("smoothx"=x,"kern_weights"=1)
        return(sx)        
    }

    
# ------------------------------------------------------------------------------


.bwNeweyWest2 <- 
function(x, kernel = c("Bartlett", "Parzen", 
    "Quadratic Spectral", "Truncated", "Tukey-Hanning"), 
    prewhite = 1, ar.method = "ols",...) 
{
    kernel <- match.arg(kernel)
    if (kernel %in% c("Truncated", "Tukey-Hanning")) 
        stop(paste("Automatic bandwidth selection only available for ", 
            dQuote("Bartlett"), ", ", dQuote("Parzen"), " and ", 
            dQuote("Quadratic Spectral"), " kernel. Use ", sQuote(".bwAndrews2"), 
            " instead.", sep = ""))
    prewhite <- as.integer(prewhite)
    n <- nrow(x)
    k <- ncol(x)
    weights <- rep(1, k)
    if (length(weights) < 2) 
        weights <- 1
    mrate <- switch(kernel, Bartlett = 2/9, Parzen = 4/25, "Quadratic Spectral" = 2/25)
    m <- floor(ifelse(prewhite > 0, 3, 4) * (n/100)^mrate)
    if (prewhite > 0) {
        x <- as.matrix(na.omit(ar(x, order.max = prewhite, 
            demean = FALSE, aic = FALSE, method = ar.method)$resid))
        n <- n - prewhite
    }
    hw <- x %*% weights
    sigmaj <- function(j) sum(hw[1:(n - j)] * hw[(j + 1):n])/n
    sigma <- sapply(0:m, sigmaj)
    s0 <- sigma[1] + 2 * sum(sigma[-1])
    s1 <- 2 * sum(1:m * sigma[-1])
    s2 <- 2 * sum((1:m)^2 * sigma[-1])
    qrate <- 1/(2 * ifelse(kernel == "Bartlett", 1, 2) + 1)
    rval <- switch(kernel, Bartlett = {
        1.1447 * ((s1/s0)^2)^qrate
    }, Parzen = {
        2.6614 * ((s2/s0)^2)^qrate
    }, "Quadratic Spectral" = {
        1.3221 * ((s2/s0)^2)^qrate
    })
    rval <- rval * (n + prewhite)^qrate
    return(rval)
}


# ------------------------------------------------------------------------------


.summary.gel <- 
function(object,interval=FALSE, ...)
{
    z <- object
    n <- nrow(z$gt)
    khat <- crossprod(z$gt)/n
    gbar <- colMeans(z$gt)
    
    se_par <- sqrt(diag(z$vcov_par))
    par <- z$par
    tval <- par/se_par

    se_parl <- sqrt(diag(z$vcov_lambda))
    lamb <- z$lambda
    tvall <- lamb/se_parl

    LR_test <- 2*z$objective*n
    LM_test <- n*crossprod(z$lambda,crossprod(khat,z$lambda))
    J_test <- n*crossprod(gbar,solve(khat,gbar))
    test <- c(LR_test,LM_test,J_test)
    vptest <- pchisq(test,(ncol(z$gt)-length(z$par)),lower.tail=FALSE)
    ans <- list(type=z$type)
    names(ans$type) <-"Type of GEL"
    
    ans$par <- round(cbind(par,se_par, tval, 2 * pnorm(abs(tval), lower.tail = FALSE)),5)
    ans$lambda <- round(cbind(lamb,se_parl, tvall, 2 * pnorm(abs(tvall), lower.tail = FALSE)),5)

        dimnames(ans$par) <- list(names(z$par), 
        c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
        dimnames(ans$lambda) <- list(names(z$lambda), 
        c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))

    ans$test <- cbind(test,vptest)
    dimnames(ans$test) <- list(c("LR test","LM test","J test"),c("statistics","p-value"))    

    if (interval != FALSE)
        {
        zs <- qnorm((1-interval)/2,lower.tail=FALSE)
        ch <- zs*se_par
        ans$interval_par <- cbind(par-ch,par+ch)
        dimnames(ans$interval_par) <- list(names(par),c("Theta_lower","Theta_upper"))
        chl <- zs*se_parl
        ans$interval_lam <- cbind(lamb-chl,lamb+chl)
        dimnames(ans$interval_lam) <- list(names(lamb),c("Lambda_lower","Lambda_upper"))
        }


    if (z$type == "EL")
        ans$badrho <- z$badrho
    if (!is.null(z$weights))
        {
        ans$weights <- z$weights
        }
    ans$conv_par <- z$conv_par
    ans$conv_pt <- z$conv_pt
    ans$conv_moment <- cbind(z$conv_moment)
    ans$conv_lambda <- z$conv_lambda
    names(ans$conv_par) <- "Convergence_code_theta"
    names(ans$conv_pt) <- "Sum_of_pt"
    names(ans$conv_lambda) <- "Convergence_code_for_lambda"
    dimnames(ans$conv_moment) <- list(names(z$gt),"Sample_moment_with_pt")
    class(ans) <- "summary.gel"
    ans    
}


# ------------------------------------------------------------------------------


.get_dat <- 
function(formula,h,intercept=TRUE) 
{
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    
    h <- cbind(rep(1,nrow(h)),h)
    colnames(h) <- c("Intercept",paste("h",1:(ncol(h)-1),sep=""))
    y <- as.matrix(model.response(mf, "numeric"))
    x <- as.matrix(model.matrix(mt, mf, NULL))
    if (!intercept)
        {
        x <- as.matrix(x[,2:ncol(x)])
        h <- as.matrix(h[,2:ncol(h)])
        }
    ny <- ncol(y)
    k <- ncol(x)
    nh <- ncol(h)
    if (nrow(y) != nrow(x) | nrow(x) != nrow(h) | nrow(y)!=nrow(h))
        stop("The number of observations of X, Y and H must be the same")
    if (nh<k)
        stop("The number of moment conditions must be at least equal to the number of coefficients to estimate")
    if (is.null(colnames(y)))
        {
        if (ny>1) 
            colnames(y) <- paste("y",1:ncol(y),sep="")
        if (ny == 1) 
            colnames(y) <- "y"
        }
    x <- cbind(y,x,h)
    return(list(x=x,nh=nh,ny=ny,k=k))
}


################################################################################

