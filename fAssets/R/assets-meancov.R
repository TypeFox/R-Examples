
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
# FUNCTION:                   DESCRIPTION:
#  assetsMeanCov               Estimates mean and variance for a set of assets
# FUNCTION:                   DESCRIPTION:
#  .covMeanCov                 uses sample covariance estimation
#  .mveMeanCov                 uses "cov.mve" from [MASS]
#  .mcdMeanCov                 uses "cov.mcd" from [MASS]
#  .studentMeanCov             uses "cov.trob" from [MASS]
#  .MCDMeanCov                 requires "covMcd" from [robustbase]  
#  .OGKMeanCov                 requires "covOGK" from [robustbase] 
#  .nnveMeanCov                uses builtin from [covRobust]
#  .shrinkMeanCov              uses builtin from [corpcor]
#  .baggedMeanCov              uses builtin from [corpcor]
#  .arwMeanCov                 uses builtin from [mvoutlier]
#  .donostahMeanCov            uses builtin from [robust]
#  .bayesSteinMeanCov          copy from Alexios Ghalanos
#  .ledoitWolfMeanCov          uses builtin from [tawny]
#  .rmtMeanCov                 uses builtin from [tawny]
# FUNCTION:                   DESCRIPTION:
#  getCenterRob                Extracts the robust estimate for the center
#  getCovRob                   Extracts the robust estimate for the covariance
################################################################################


assetsMeanCov <- 
  function(x, 
           method = c("cov", "mve", "mcd", "MCD", "OGK", "nnve", "shrink", "bagged"), 
           check = TRUE, force = TRUE, baggedR = 100, sigmamu = scaleTau2, alpha = 1/2,
           ...)
  {   
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Computes robust mean and covariance from multivariate time series
    
    # Arguments:
    #   x - a multivariate time series, a data frame, or any other
    #       rectangular object of assets which can be converted into
    #       a matrix by the function 'as.matrix'. Optional Dates are 
    #       rownames, instrument names are column names.
    #   method - Which method should be used to compute the covarinace?
    #       method = "cov"        sample covariance computation
    #       method = "mve"        uses "mve" from [MASS]
    #       method = "mcd"        uses "mcd" from [MASS]
    #       method = "MCD"        uses "MCD" from [robustbase]
    #       method = "OGK"        uses "OGK" from [robustbase]
    #       method = "nnve"       uses "nnve" from [covRobust]
    #       method = "shrink"     uses "shrinkage" from [corpcor] 
    #       method = "bagged"     uses "bagging" [corpcor]
    #   alpha - MCD: numeric parameter controlling the size of the subsets 
    #       over which the determinant is minimized, i.e., alpha*n observations 
    #       are used for computing the determinant. Allowed values are between 
    #       0.5 and 1 and the default is 0.5.
    #   sigma.mu - OGK: a function that computes univariate robust location 
    #       and scale estimates. By default it should return a single numeric 
    #       value containing the robust scale (standard deviation) estimate. 
    #       When mu.too is true, sigmamu() should return a numeric vector of 
    #       length 2 containing robust location and scale estimates. See 
    #       scaleTau2, s_Qn, s_Sn, s_mad or s_IQR for examples to be used as 
    #       sigmamu argument.
    
    # Note:
    #   The output of this function can be used for portfolio
    #   optimization.
    
    # Example:
    #   DJ = 100 * returns(as.timeSeries(data(DowJones30)))
    #   DJ = DJ[, c("CAT", "IBM", "GE", "JPM")]
    #   Sample Covariance:
    #       assetsMeanCov(DJ, "cov")
    #   MASS:
    #       assetsMeanCov(DJ, "mve")
    #       assetsMeanCov(DJ, "mcd")
    #   require(robustbase)
    #       assetsMeanCov(DJ, "MCD")
    #       assetsMeanCov(DJ, "OGK")
    #   require(covRobust)
    #       assetsMeanCov(DJ, "nnve")
    
    # FUNCTION:
    
    # Transform Input:
    x.mat <- as.matrix(x)
    
    # Do not use: method = match.arg(method)
    method <- method[1]
    N <- ncol(x)
    assetNames <- colnames(x.mat)
    
    # Attribute Control List:
    control <- c(method = method[1])
    user <- TRUE
    
    # Compute Classical Covariance:
    if (method == "cov") {
      # Classical Covariance Estimation:
      ans = list(center = colMeans(x.mat), cov = cov(x.mat))
      user = FALSE
    }
    
    # From R Package "robustbase":
    if (method == "MCD" | method == "Mcd") {
      ans <- robustbase::covMcd(x.mat, alpha = alpha, ...)
      mu = ans$center
      Sigma = ans$cov
      user = FALSE
    }   
    if (method == "OGK" | method == "Ogk") {
      ans <- robustbase::covOGK(x.mat, sigmamu = sigmamu, ...)
      user = FALSE
    }
    
    # [MASS] mve and mcd Routines:
    if (method == "mve") {
      ans = MASS::cov.rob(x = x.mat, method = "mve", ...)
      user = FALSE
    }
    if (method == "mcd") {
      ans = MASS::cov.rob(x = x.mat, method = "mcd", ...) 
      user = FALSE
    }    
    
    # [corpcor] Shrinkage and Bagging Routines 
    if (method == "shrink") {
      fit = .cov.shrink(x = x.mat, ...)
      ans = list(center = colMeans(x.mat), cov = fit)
      user = FALSE
    } 
    if (method == "bagged") {
      fit = .cov.bagged(x = x.mat, R = baggedR, ...)
      ans = list(center = colMeans(x.mat), cov = fit)
      control = c(control, R = as.character(baggedR))
      user = FALSE
    }
    
    # Nearest Neighbour Variance Estimation:
    if (method == "nnve") {
      fit = .cov.nnve(datamat = x.mat, ...)
      ans = list(center = colMeans(x.mat), cov = fit$cov)
      user = FALSE
    }
    
    # User specified estimator:
    if(user) {
      fun = match.fun(method[1])
      ans = fun(x.mat, ...)
    }
    
    # Result:
    mu = center = ans$center
    Sigma = cov = ans$cov
    
    # Add Size to Control List:
    control = c(control, size = as.character(N))
    
    # Add Names for Covariance Matrix to Control List:
    names(mu) = assetNames
    colnames(Sigma) = rownames(Sigma) = colNames = assetNames
    
    # Check Positive Definiteness:
    if (check) {
      result = isPositiveDefinite(Sigma)
      if(result) {
        control = c(control, posdef = "TRUE")
      } else {
        control = c(control, posdef = "FALSE")
      }
    }
    
    # Check Positive Definiteness:
    control = c(control, forced = "FALSE")
    if (force) {
      control = c(control, forced = "TRUE")
      if (!result) Sigma = makePositiveDefinite(Sigma)       
    }
    
    # Result:
    ans = list(center = mu, cov = Sigma, mu = mu, Sigma = Sigma)
    attr(ans, "control") = control
    
    # Return Value:
    ans
  }


################################################################################


.covMeanCov <- 
  function(x, ...)
  {
    # Description:
    #   Uses sample covariance estimation
    
    # Arguments:
    #   x - an object of class timeSeries
    
    # FUNCTION:
    
    # Settings:
    x.mat = as.matrix(x)
    N = ncol(x)
    assetNames = colnames(x)
    
    ans = list(center = colMeans(x.mat), cov = cov(x.mat))
    names(ans$center) = assetNames
    rownames(ans$cov) = colnames(ans$cov) = assetNames 
    
    # Return Value:
    ans
  }


# ------------------------------------------------------------------------------


.mveMeanCov <- 
  function(x, ...)
  {
    # Description:
    
    # Arguments:
    #   x - an object of class timeSeries
    
    # FUNCTION:
    
    # Settings:
    x.mat = as.matrix(x)
    N = ncol(x)
    assetNames = colnames(x)
    
    ans <- MASS::cov.rob(x = x.mat, method = "mve")
    names(ans$center) = assetNames
    rownames(ans$cov) = colnames(ans$cov) = assetNames
    
    # Return Value:
    ans
  }


# ------------------------------------------------------------------------------


.mcdMeanCov <- 
  function(x, ...)
  {
    # Description:
    
    # Arguments:
    #   x - an object of class timeSeries
    
    # FUNCTION:
    
    # Settings:
    x.mat = as.matrix(x)
    N = ncol(x)
    assetNames = colnames(x)
    
    ans <- MASS::cov.rob(x = x.mat, method = "mcd") 
    names(ans$center) = assetNames
    rownames(ans$cov) = colnames(ans$cov) = assetNames
    
    # Return Value:
    ans
  }


# ------------------------------------------------------------------------------


.studentMeanCov <-
  function(x, ...)
  {
    # Description:
    
    # Arguments:
    #   x - an object of class timeSeries
    
    # FUNCTION:
    
    # Settings:
    x.mat = as.matrix(x)
    N = ncol(x)
    assetNames = colnames(x)
    
    ans <- MASS::cov.trob(x, ...)
    names(ans$center) = assetNames
    rownames(ans$cov) = colnames(ans$cov) = assetNames
    
    # Return Value:  
    ans
  }


# ------------------------------------------------------------------------------


.MCDMeanCov <- 
  function(x, alpha = 1/2, ...)
  {
    # Description:
    
    # Arguments:
    #   x - an object of class timeSeries
    
    # FUNCTION:
    
    # Settings:
    x.mat = as.matrix(x)
    N = ncol(x)
    assetNames = colnames(x)
    
    ans <- robustbase::covMcd(x.mat, alpha = alpha, ...)
    names(ans$center) = assetNames
    rownames(ans$cov) = colnames(ans$cov) = assetNames
    
    # Return Value:
    ans
  }


# ------------------------------------------------------------------------------


.OGKMeanCov <- 
  function(x, sigmamu = scaleTau2, ...)
  {
    # Description:
    
    # Arguments:
    #   x - an object of class timeSeries
    
    # FUNCTION:
    
    # Settings:
    x.mat = as.matrix(x)
    N = ncol(x)
    assetNames = colnames(x)
    
    ans <- robustbase::covOGK(x.mat, sigmamu = sigmamu, ...)
    names(ans$center) = assetNames
    rownames(ans$cov) = colnames(ans$cov) = assetNames
    
    # Return Value:
    ans    
  }


# ------------------------------------------------------------------------------


.nnveMeanCov <- 
  function(x, ...)
  {
    # Description:
    
    # Arguments:
    #   x - an object of class timeSeries
    
    # FUNCTION:
    
    # Settings:
    x.mat = as.matrix(x)
    N = ncol(x)
    assetNames = colnames(x)
    fit = .cov.nnve(datamat = x.mat, ...)
    
    ans = list(center = colMeans(x.mat), cov = fit$cov) 
    names(ans$center) = assetNames
    rownames(ans$cov) = colnames(ans$cov) = assetNames
    
    # Return Value:
    ans
  }


# ------------------------------------------------------------------------------


.shrinkMeanCov <- 
  function(x, ...)
  {
    # Description:
    
    # Arguments:
    #   x - an object of class timeSeries
    
    # Note:                                              
    #   Based on a function borrowed from package corpcor
    
    # FUNCTION:
    
    # Settings:
    x.mat = as.matrix(x)
    N = ncol(x)
    assetNames = colnames(x)
    fit = .cov.shrink(x = x.mat, ...)
    
    ans = list(center = colMeans(x.mat), cov = fit)
    names(ans$center) = assetNames
    rownames(ans$cov) = colnames(ans$cov) = assetNames
    
    # Return Value:
    ans
  }


# ------------------------------------------------------------------------------


.baggedMeanCov <- 
  function(x, baggedR = 100, ...)
  {
    # Description:
    
    # Arguments:
    #   x - an object of class timeSeries
    
    # Note:                                                
    #   Based on a function borrowed from package corpcor
    
    # FUNCTION:
    
    # Settings:
    x.mat = as.matrix(x)
    N = ncol(x)
    assetNames = colnames(x)
    fit = .cov.bagged(x = x.mat, R = baggedR, ...)
    
    ans = list(center = colMeans(x.mat), cov = fit)
    names(ans$center) = assetNames
    rownames(ans$cov) = colnames(ans$cov) = assetNames
    
    # Return Value:
    ans
  }


# ------------------------------------------------------------------------------


.arwMeanCov <- 
  function(x, ...)
  {
    # Description:
    #   Adaptive reweighted estimator for multivariate location and scatter
    #   with hard-rejection weights and delta = chi2inv(1-d,p)
    
    # Arguments:
    #   x - an object of class timeSeries
    
    # Note:
    #   Based on a function borrowed from package mvoutlier
    
    # FUNCTION:
    
    # Settings:
    x.mat = as.matrix(x)
    N = ncol(x)
    assetNames = colnames(x)
    fit = .cov.arw(x = x.mat, center = colMeans(x.mat), cov = cov(x),, ...)
    
    ans = list(center = fit$center, cov = fit$cov)
    names(ans$center) = assetNames
    rownames(ans$cov) = colnames(ans$cov) = assetNames
    
    # Return Value:
    ans
  }


# ------------------------------------------------------------------------------


.donostahMeanCov <- 
  function(x, ...)
  {
    # Description:
    
    # Arguments:
    #   x - an object of class timeSeries
    
    # Note:
    #   Based on a function borrowed from package robust
    
    # Settings:
    x.mat = as.matrix(x)
    N = ncol(x)
    assetNames = colnames(x)
    
    ans = .cov.donostah(x = x.mat)
    names(ans$center) = assetNames
    rownames(ans$cov) = colnames(ans$cov) = assetNames
    
    # Return Value:
    ans
  }


# ------------------------------------------------------------------------------


.bayesSteinMeanCov <- 
  function(x, ...)
  {
    # Description:
    
    # Arguments:
    #   x - an object of class timeSeries
    
    # Note:
    #   Based on a function written by Alexios Ghalanos
    
    # Bayes Stein estimator
    # Alexios Ghalanos 2008
    # alexios at 4dscape.com
    # This function encapsulates an example of shrinking the returns 
    #   and covariance using Bayes-Stein shrinkage as described in 
    #   Jorion, 1986.
    
    # Settings:
    data <- getDataPart(x)
    mu <- as.matrix(apply(data,2, FUN = function(x) mean(x)))
    S <- cov(data)
    k <- dim(data)[2]
    n <- dim(data)[1]
    one <- as.matrix(rep(1, k))
    a <- solve(S, one)
    
    # Constant non informative prior
    mu.prior <- one * as.numeric(t(mu) %*% a/t(one) %*% a)
    S.inv <- solve(S)
    d <- t(mu-mu.prior) %*% S.inv %*% (mu-mu.prior)
    d <- as.numeric(d)
    lambda <- (k+2) / d
    w <- lambda / (n+lambda)
    mu.pred <- (1-w) * mu + w * mu.prior
    
    wc1 <- 1 / (n+lambda)
    wc2 <- lambda*(n-1) / (n*(n+lambda)*(n-k-2))
    wc2 <- wc2 / as.numeric(t(one) %*% a)
    V.post <- wc1 * S + wc2 * one %*% t(one)
    V.pred <- S + V.post
    sigma.post <- sqrt(diag(V.post))
    sigma.pred <- sqrt(diag(V.pred))
    
    result <- list(
      mu = mu, mu.prior = mu.prior, mu.predict = mu.pred, 
      V = S, V.post = V.post, V.pred = V.pred, Sigma = sqrt(diag(S)), 
      Sigma.post = sigma.post, Sigma.predict = sigma.pred)
    
    ans = list(center = result$mu.pred[,1], cov = result$V.pred)
    names(ans$center) = colnames(x)
    rownames(ans$cov) = colnames(ans$cov) = colnames(x)
    
    # Return Value:
    ans
  }


# ------------------------------------------------------------------------------


.ledoitWolfMeanCov <-
  function(x, ...)
  {
    # Description:
    #   Perform shrinkage on a sample covariance towards a biased covariance
    
    # Arguments:
    #   x - an object of class timeSeries
    
    # Details:
    #   This performs a covariance shrinkage estimation as specified in 
    #   Ledoit and Wolf. Using within the larger framework only requires 
    #   using the getCorFilter.Shrinkage function, which handles the work 
    #   of constructing a shrinkage estimate of the covariance matrix of 
    #   returns (and consequently its corresponding correlation matrix).
    
    # Note:
    #   Based on a function borrowed from package tawny
    #   Author: Brian Lee Yung Rowe
    
    # Settings:
    data = getDataPart(x)
    center = colMeans(data)
    cov = .cov.shrink.tawny(data, ...)
    
    ans = list(center = center, cov = cov)
    names(ans$center) = colnames(x)
    rownames(ans$cov) = colnames(ans$cov) = colnames(x)
    
    # Return Value:
    ans
  }


# ------------------------------------------------------------------------------


.rmtMeanCov <-
  function(x, ...)
  {
    # Description:
    #   Perform Random Matrix Theory on correlation matrix
    
    # Arguments:
    #   x - an object of class timeSeries
    
    # Author: 
    #   tawnyBrian Lee Yung Rowe
    
    # Note:
    #   Based on a function borrowed from package tawny
    #   Author: Brian Lee Yung Rowe
    
    # FUNCTION:
    
    # Settings:
    data = getDataPart(x)
    center = colMeans(data)
    cor = .filter.RMT(data, trace = FALSE, doplot = FALSE)
    
    g = colSds(data)
    N = length(g)
    cov = 0*cor
    for (i in 1:N)
      for (j in i:N)
        cov[i,j] = cov[j,i] = g[i]*cor[i,j]*g[j]
    
    ans = list(center = center, cov = cov)
    names(ans$center) = colnames(x)
    rownames(ans$cov) = colnames(ans$cov) = colnames(x)
    
    # Return Value:
    ans
  }


################################################################################


getCenterRob <-
  function(object)
  {
    # Return Value:
    object$center
  }


# ------------------------------------------------------------------------------


getCovRob <-
  function(object)
  {
    # Return Value:
    object$cov
  }


################################################################################

