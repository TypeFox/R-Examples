# fields  is a package for analysis of spatial data written for
# the R software environment .
# Copyright (C) 2016
# University Corporation for Atmospheric Research (UCAR)
# Contact: Douglas Nychka, nychka@ucar.edu,
# National Center for Atmospheric Research, PO Box 3000, Boulder, CO 80307-3000
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the R software environment if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
# or see http://www.r-project.org/Licenses/GPL-2    
MaternGLS.test <- function(x, y, smoothness = 1.5, 
    init = log(c(1, 0.2, 0.1))) {
    # some simulations within fields to
    # study the variability in estimates of the covariance parameters.
    N <- length(y)
    Tmatrix <- fields.mkpoly(x, 2)
    qr.T <- qr(Tmatrix)
    Q2 <- qr.yq2(qr.T, diag(1, N))
    nu <- smoothness
    
    loglmvn <- function(pars, nu, y, x, d) {
        lrho = pars[1]
        ltheta = pars[2]
        lsig2 = pars[3]
        #  print( pars)
        N <- length(y)
        M <- (exp(lrho) * Matern(d, range = exp(ltheta), smoothness = nu) + 
            exp(lsig2) * diag(N))
        
        X <- fields.mkpoly(x, 2)
        Mi <- solve(M)
        betahat <- solve(t(X) %*% Mi %*% X) %*% t(X) %*% Mi %*% 
            y
        res <- y - X %*% betahat
        cM <- chol(M)
        lLike <- (N/2) * log(2 * pi) - (1/2) * (2 * sum(log(diag(cM)))) - 
            (1/2) * t(res) %*% Mi %*% res
        
        ycept <- -lLike
        
        if ((abs(lrho) > 20) | (abs(ltheta) > 10) | (abs(lsig2) > 
            40)) {
            return(ycept + 1000 * sum(abs(abs(pars) - 100)))
        }
        else {
            return(ycept)
        }
    }
    
    d <- rdist(x, x)
    temp <- optim(init, loglmvn, method = "L-BFGS-B", nu = nu, 
        y = y, x = x, d = d)
    out <- exp(temp$par)
    
    return(list(smoothness = smoothness, pars = out, optim = temp))
}



MaternGLSProfile.test <- function(x, y, smoothness = 1.5, 
    init = log(c(0.05, 1))) {
    # some simulations within fields to
    # study the variability in estimates of the covariance parameters.
    N <- length(y)
    
    nu <- smoothness
    
    loglmvn <- function(pars, nu, y, x, d) {
        llam = pars[1]
        ltheta = pars[2]
        # print( pars)
        N <- length(y)
        theta <- exp(ltheta)
        lambda <- exp(llam)
        lLike <- mKrig(x, y, theta = theta, Covariance = "Matern", 
            smoothness = nu, lambda = lambda)$lnProfileLike
        ycept <- -lLike
        if ((abs(llam) > 20) | (abs(ltheta) > 10)) {
            return(ycept + 1000 * sum(abs(abs(pars) - 100)))
        }
        else {
            return(ycept)
        }
    }
    
    d <- rdist(x, x)
    temp <- optim(init, loglmvn, method = "L-BFGS-B", nu = nu, 
        y = y, x = x, d = d)
    out <- exp(temp$par)
    rho.MLE <- mKrig(x, y, theta = out[2], Covariance = "Matern", 
        smoothness = nu, lambda = out[1])$rho.MLE
    sigma2.MLE <- out[1] * rho.MLE
    return(list(smoothness = smoothness, pars = c(rho.MLE, out[2], 
        sigma2.MLE), optim = temp))
}

MaternQR.test <- function(x, y, smoothness = 1.5, 
    init = log(c(1, 0.2, 0.1))) {
    # some simulations within fields to
    # study the variability in estimates of the covariance parameters.
    nu <- smoothness
    
    loglmvn <- function(pars, nu, x, y) {
        N <- length(y)
        Tmatrix <- fields.mkpoly(x, 2)
        qr.T <- qr(Tmatrix)
        Q2 <- qr.yq2(qr.T, diag(1, N))
        ys <- t(Q2) %*% y
        N2 <- length(ys)
        lrho = pars[1]
        ltheta = pars[2]
        lsig2 = pars[3]
        d <- rdist(x, x)
        
        A <- (exp(lrho)*Matern(d, range = exp(ltheta), 
            smoothness = nu) + exp(lsig2) * diag(N))
        A <- t(Q2) %*% A %*% Q2
        A <- chol(A)
        w = backsolve(A, ys, transpose = TRUE)
        ycept <- (N2/2) * log(2 * pi) + sum(log(diag(A))) + (1/2) * 
            t(w) %*% w
        
        if ((abs(lrho) > 100) | (abs(ltheta) > 100) | (abs(ltheta) > 
            100)) {
            return(ycept + 1000 * sum(abs(abs(pars) - 100)))
        }
        else {
            return(ycept)
        }
    }
    
    
    temp <- optim(init, loglmvn, method = "L-BFGS-B", nu = nu, 
        x = x, y = y)
    out <- exp(temp$par)
    llike <- loglmvn(temp$par, nu, x, y)
    return(list(smoothness = smoothness, pars = out, llike = llike, 
        optim = temp))
}


MaternQRProfile.test <- function(x, y, smoothness = 1.5, 
    init = log(c(1))) {
    # some simulations within fields to
    # study the variability in estimates of the covariance parameters.
    nu <- smoothness
    loglmvn <- function(pars, nu, x, y) {
        ltheta = pars[1]
        # print( exp(ltheta))
        ycept <- Krig(x, y, Covariance = "Matern", theta = exp(ltheta), 
            smoothness = nu, method = "REML")$lambda.est[6, 5]
        #        print( c(exp(ltheta),ycept))
        if ((abs(ltheta) > 100)) {
            return(ycept + 1000 * sum(abs(abs(pars) - 100)))
        }
        else {
            return(ycept)
        }
    }
    
    temp <- optim(init, loglmvn, method = "L-BFGS-B", nu = nu, 
        x = x, y = y)
    theta.est <- exp(temp$par[1])
    out2 <- Krig(x, y, Covariance = "Matern", theta = theta.est, 
        smoothness = nu, method = "REML")
    # MLE based on reduced degrees of freedom:
    
    offset <- (out2$N/(out2$N - 3))
    
    out3 <- c(out2$rho.MLE * offset, theta.est, out2$shat.MLE^2 * 
        offset)
    
    return(list(obj = out2, smoothness = smoothness, pars = out3, 
        trA = out2$eff.df, optim = temp))
}

# this function has correct formula for REML likelihood
REML.test <- function(x, y, rho, sigma2, theta, nu = 1.5) {
    Tmatrix <- fields.mkpoly(x, 2)
    qr.T <- qr(Tmatrix)
    N <- length(y)
    Q2 <- qr.yq2(qr.T, diag(1, N))
    ys <- t(Q2) %*% y
    N2 <- length(ys)
    A <- (rho * Matern(rdist(x, x), range = theta, smoothness = nu) + 
        sigma2 * diag(1, N))
    A <- t(Q2) %*% A %*% Q2
    Ac <- chol(A)
    w <- backsolve(Ac, ys, transpose = TRUE)
    REML.like <- (N2/2) * log(2 * pi) + (1/2) * 2 * sum(log(diag(Ac))) + 
        (1/2) * t(w) %*% w
    REML.like <- -1 * REML.like
    ccoef <- rho * Q2 %*% solve(A) %*% ys
    return(list(REML.like = REML.like, A = A, ccoef = ccoef, 
        quad.form = t(w) %*% w, rhohat = (t(w) %*% w/N2) * rho, 
        det = 2 * sum(log(diag(Ac))), N2 = N2))
}


