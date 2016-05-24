####**********************************************************************
####**********************************************************************
####
####  BOOSTED MULTIVARIATE TREES FOR LONGITUDINAL DATA (BOOSTMTREE)
####  Version 1.1.0 (_PROJECT_BUILD_ID_)
####
####  Copyright 2016, University of Miami
####
####  This program is free software; you can redistribute it and/or
####  modify it under the terms of the GNU General Public License
####  as published by the Free Software Foundation; either version 3
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
####  Project Partially Funded By:
####  ----------------------------------------------------------------
####  Dr. Ishwaran's work was funded in part by grant R01 CA163739 from
####  the National Cancer Institute.
####
####  Dr. Kogalur's work was funded in part by grant R01 CA163739 from 
####  the National Cancer Institute.
####  ----------------------------------------------------------------
####  Written by:
####  ----------------------------------------------------------------
####    Hemant Ishwaran, Ph.D.
####    Professor, Division of Biostatistics
####    Clinical Research Building, Room 1058
####    1120 NW 14th Street
####    University of Miami, Miami FL 33136
####
####    email:  hemant.ishwaran@gmail.com
####    URL:    http://web.ccs.miami.edu/~hishwaran
####    --------------------------------------------------------------
####    Amol Pande
####    Division of Biostatistics
####    1120 NW 14th Street
####    University of Miami, Miami FL 33136
####
####    email:  amoljpande@gmail.com
####    --------------------------------------------------------------
####    Udaya B. Kogalur, Ph.D.
####    Consultant Staff
####    Deptartment of Quantitative Health Sciences
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


simLong <- function(n = 100,
                     ntest = 0,
                     N = 5,
                     rho = 0.8,
                     type = c("corCompSym", "corAR1", "corSymm", "iid"),
                     model = c(0, 1, 2, 3),
                     phi = 1,
                     q = 0,
                     ...)
{
  type <- match.arg(type, c("corCompSym", "corAR1", "corSymm", "iid"))
  model <- as.numeric(match.arg(as.character(model), as.character(0:3)))
  dta <- data.frame(do.call("rbind", lapply(1:(n+ntest), function(i) {
    Ni <- round(runif(1, 1, 3 * N))
    type <- match.arg(type, c("corCompSym", "corAR1", "corSymm", "iid"))
    if (type == "corCompSym") {
      corr.mat <- matrix(rho, nrow=Ni, ncol=Ni)
      diag(corr.mat) <- 1
    }
    if (type == "corAR1") {
      corr.mat <- diag(rep(1, Ni))
      if (Ni > 1) {
       for (ii in 1:(Ni - 1)) {
          corr.mat[ii, (ii + 1):Ni] <- rho^(1:(Ni - ii))
        }
        ind <- lower.tri(corr.mat) 
        corr.mat[ind] <- t(corr.mat)[ind]
      }
    }
    if (type == "iid") {
      corr.mat <- diag(rep(1, Ni))
    }
    eps <- sqrt(phi) * t(chol(corr.mat)) %*% rnorm(Ni)
    x1 <- rnorm(1)
    x2 <- runif(1, 1, 2)
    x3 <- runif(1, 1, 3)
    x4 <- rnorm(1)
    x <- c(x1, x2, x3, x4)
    p <- length(x)
    if (q > 0) {
      xnoise <- rnorm(q)
      x <- c(x, xnoise)
    }
    tm <- sample((1:(3 * N))/N, size = Ni, replace = TRUE)
    if (model == 0) {
      y <- 1.5 + 2.5 * x1 - 1.2 * x3 - .6 * x4 + eps
    }
    if (model == 1) {
      y <- 1.5 + 2.5 * x1 - 1.2 * x3 - .2 * x4 - .65 * tm  * x2   + eps
    }
    if (model == 2) {
      y <- 1.5 + 2.5 * x1 - 1.2 * x3 - .2 * x4 - .65 * (tm ^ 2) * (x2 ^ 2)   + eps
    }
    if (model == 3) {
      y <- 1.5 + 2.5 * x1 - 1.2 * x3 - .2 * exp(x4) - .65 * (tm ^ 2) * (x2 ^2) * x3  + eps
    }
    cbind(matrix(x, nrow = Ni, ncol = length(x), byrow = TRUE),
          tm, rep(i, Ni), y)
  })))
  d <- q + 4
  colnames(dta) <- c(paste("x", 1:d, sep = ""), "time", "id", "y")
  dtaL <- list(features = dta[, 1:d], time = dta$time, id = dta$id, y = dta$y) 
  if (model == 0) {
    f.true <- "y ~ x1 + x3 + x4"
  }
  if (model == 1) {
    f.true <- "y ~ x1 + x3 + x4 + I(time * x2)"
  }
  if (model == 2) {
    f.true <- "y ~ x1 + x3 + x4 + I(time^2 * x2^2)"
  }
  if (model == 3) {
    f.true <- "y ~ x1 + x3 + exp(x4) + I(time^2 * x2^2 * x3)"
  }
  trn <- c(1:sum(dta$id <= n))
  return(invisible(list(dtaL = dtaL, dta = dta, trn = trn, f.true = f.true)))
}
