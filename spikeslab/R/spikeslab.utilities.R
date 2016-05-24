####**********************************************************************
####**********************************************************************
####
####  SPIKE AND SLAB 1.1.5
####
####  Copyright 2013, Cleveland Clinic Foundation
####
####  This program is free software; you can redistribute it and/or
####  modify it under the terms of the GNU General Public License
####  as published by the Free Software Foundation; either version 2
####  of the License, or (at your option) any later version.
####
####  This program is distributed in the hope that it will be useful,
####  but WITHOUT ANY WARRANTY; without even the implied warranty of
####  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
####  GNU General Public License for more details.
####
####  You should have received a copy of the GNU General Public
####  License along with this program; if not, write to the Free
####  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
####  Boston, MA  02110-1301, USA.
####
####  ----------------------------------------------------------------
####  Written and Developed by:
####  ----------------------------------------------------------------
####    Hemant Ishwaran, Ph.D.
####    Director of Statistical Methodology
####    Professor, Division of Biostatistics
####    Clinical Research Building, Room 1058
####    1120 NW 14th Street
####    University of Miami, Miami FL 33136
####
####    email:  hemant.ishwaran@gmail.com
####    URL:    http://web.ccs.miami.edu/~hishwaran
####    --------------------------------------------------------------
####    J. Sunil Rao, Ph. D.
####    Professor and Director of the Division of Biostatistics, 
####    Department of Epidemiology & Public Health
####    Clinical Research Bldg, R-669
####    1120 NW 14th Street, Room 1056
####    Miami, FL 33136
####    email:  rao.jsunil@gmail.com
####    URL:    http://biostat.med.miami.edu/people/primary-faculty/sunil-rao
####  ----------------------------------------------------------------
####  Maintained by:
####    Udaya B. Kogalur, Ph.D.
####    Adjunct Staff
####    Dept of Quantitative Health Sciences
####    Cleveland Clinic Foundation
####    
####    Kogalur & Company, Inc.
####    5425 Nestleway Drive, Suite L1
####    Clemmons, NC 27012
####
####    email:  ubk@kogalur.com
####    URL:    http://www.kogalur.com
####    --------------------------------------------------------------
####
####**********************************************************************
####**********************************************************************

### --------------------------------------------------------------
###
###  spike slab core internal functions
###
###--------------------------------------------------------------


# nice standard errors for plots (shamelessly taken from lars)
error.bars <- function(x, upper, lower, width = 0.0025, ...)
{
        xlim <- range(x)
        barw <- diff(xlim) * width
        segments(x, upper, x, lower, ...)
        segments(x - barw, upper, x + barw, upper, ...)
        segments(x - barw, lower, x + barw, lower, ...)
}

# sd for rescaling X
sd.center <- function(x, center) {
  n.x <- length(x)
  if (center) {
    sd(x, na.rm=T)*sqrt((n.x-1)/n.x)
  }
  else {
    sqrt(mean(x^2, na.rm=T))
  }
}

# mean for centering X
mean.center <- function(x, center) {
  if (center) {
    mean(x, na.rm=T)
  }
  else {
    0
  }
}

# robust SD
SD <- function(x) {
  if (all(is.na(x))) return(NA) else return(sd(x, na.rm=T))
}

# robust MSE
MSE <- function(x) {
  median(x^2, na.rm=T)
}

# XX matrix
XX.multiply <- function(A, bigp.smalln) {
  if (bigp.smalln) {
    NULL 
  }
  else {
    as.matrix(t(A)%*%A)
  }   
}

# X %*% b calculation when many b's are zero
X.b.mult <- function(A, b, pt) {
  as.matrix(A[, pt]) %*% b[pt]
}

# svdwrapper
# work around for La.svd(x, nu, nv): error code 1 from Lapack routine 'dgesdd'
# tries svd(t(x)) in case of failure: suggested by Art Owen (R-help)
svdwrapper <- function(x) {
  gotit <- F
  try({svdx = svd(x); gotit = TRUE}, silent = TRUE)
  if(gotit) return(svdx)
  try({svdtx = svd(t(x)); gotit = TRUE}, silent = TRUE)
  if(!gotit) stop("svd(x) and svd(t(x)) both failed")
  warning("svd(x) failed but svd(t(x)) worked")
  temp <- svdtx$u
  svdtx$u <- svdtx$v
  svdtx$v <- temp
  svdtx
}

# nice wrapper for ridge
Ridge <- function(Y, X, l) {
    n <- nrow(X)
    p <- ncol(X)
    Xs <- svdwrapper(X)
    rhs <- t(Xs$u) %*% Y
    d <- Xs$d
    div <- d^2 + rep(l, length(d))
    a <- drop(d * rhs)/div
    c(Xs$v %*% a)
}

# generalized elastic net (for variable selection)
gnet.get <- function(Y, X, X.pt, phat, penal, mse = NULL, eps = .Machine$double.eps) {
  # gnet's solution uses ellipsoid optimization centered at the GRR estimator
  # most closely aligned with the BMA
  # gnet's solution set determined using AIC 
  p <- ncol(X)
  n <- nrow(X)
  if (phat == 0) {
    return(list(final.model = rep(0, p), coef = NULL))
  }
  X <- as.matrix(X[, X.pt[1:phat]])
  penal <- penal[X.pt][1:phat]
  penal[is.nan(penal)] <- Inf
  if (all(is.infinite(penal)) | all(is.na(penal))) {
    return(list(gnet = rep(0, p), gnet.path = NULL))
  }
  penal[penal < 1e-10] <- 1e-10
  penal[is.na(penal) | is.infinite(penal)] <- 1e10
  X.org <- X
  Y.org <- Y
  X <- rbind(X.org, diag(penal^0.5, ncol(X.org)))
  Y <- c(Y.org, rep(0, phat))
  lasso.out <- lars(X, Y, intercept = FALSE, use.Gram = FALSE, type = "lar")
  lasso.pred <- predict.lars(lasso.out, X.org, type = "fit", mode = "fraction")$fit
  lasso.coef <- as.matrix(predict.lars(lasso.out, type = "coef", mode = "fraction")$coef)
  dof <- apply(lasso.coef, 1, function(x) {sum(abs(x) > eps)})
  aic.penal <- sapply(1:ncol(lasso.pred), function(j) {
      mean((Y.org - lasso.pred[, j])^2) + 2 * mse * dof[j] /n
  })
  gnet <- rep(0, p)
  gnet[X.pt[1:phat]] <- lasso.coef[min(which(aic.penal == min(aic.penal, na.rm = TRUE))), ]
  gnet.path <- matrix(0, nrow(lasso.coef), p)
  gnet.path[, X.pt[1:phat]] <- lasso.coef
  gnet.path <- list(path = gnet.path, aic = aic.penal)
  return(list(gnet = gnet, gnet.path = gnet.path, gnet.object = lasso.out))
}

### --------------------------------------------------------------
###  Gibbs core internal functions
###
###  o  robust correlation
###  o  resampling (corrected for n = 1)
###  o  split for beta-blocking in Gibbs
###  o  two-component update
### --------------------------------------------------------------

Cor <- function(y,x) {
    if ((sd(y)==0) | (sd(x)==0)) {
      0
    }
    else {
      cor(y, x, use="compl")
    }  
}
  
resample <- function(x, size, ...) {
    if (length(x) <= 1) {
      if (!missing(size) && size == 0) x[FALSE] else x
    }
    else
      sample(x, size, ...)
}

beta.fold <- function(p, blocks = 5) {
    split(resample(1:p), rep(1:blocks, length = p))
}


varianceSample <- function(b,w,v0,V) {
   w1 <- (1-w)*exp(-b^2/(2*v0))/sqrt(v0)
   w2 <- w*exp(-b^2/(2*V))/sqrt(V)
   if (w1==0) w1 <- 1e-4  #numerical issue
   if (w2==0) w2 <- 1e-4  #numerical issue
   resample(c(v0,V), 1, prob=c(w1,w2))
}


 
