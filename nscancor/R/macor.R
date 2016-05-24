#  Copyright 2013 Christian Sigg
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

#' Multi-Domain Additional Explained Correlation
#' 
#' \code{macor} generalizes \code{\link{acor}} to the case of more than two data
#' domains.
#' 
#' @export
#' @param x a list of numeric matrices which contain the data from the different
#'   domains
#' @param coef a list of matrices containing the canonical vectors related to
#'   each data domain. Each matrix contains the respective canonical vectors as
#'   its columns.
#' @param center a list of logical values indicating whether the empirical mean 
#'   of (each column of) the corresponding data matrix should be subtracted. 
#'   Alternatively, a list of vectors can be supplied, where each vector 
#'   specifies the mean to be subtracted from the corresponding data matrix. 
#'   Each list element is passed to \code{\link{scale}}.
#' @param scale_ a list of logical values indicating whether the columns of the 
#'   corresponding data matrix should be scaled to have unit variance before the
#'   analysis takes place. The default is \code{FALSE} for consistency with 
#'   \code{acor}. Alternatively, a list of vectors can be supplied, where each 
#'   vector specifies the standard deviations used to rescale the columns of the
#'   corresponding data matrix. Each list element is passed to 
#'   \code{\link{scale}}.
#' @return \code{macor} returns a list of class \code{mcancor} with the 
#'   following elements: \item{cor}{a multi-dimensional array containing the 
#'   additional correlations explained by each pair of canonical variables. The 
#'   first two dimensions correspond to the domains, and the third dimension 
#'   corresponds to the different canonical variables per domain.} 
#'   \item{coef}{copied from the input arguments} \item{center}{the list of
#'   empirical means used to center the data matrices} \item{scale}{the list of
#'   empirical standard deviations used to scale the data matrices}\item{xp}{the
#'   list of deflated data matrices corresponding to \code{x}}
macor <- function(x, coef, center = TRUE, scale_ = FALSE) {
  
  X <- x
  W <- coef
  m <- length(X)  # number of domains
  n <- nrow(X[[1]])  # number of observations
  nvar <- ncol(W[[1]])  # number of canonical variables for each domain
  
  Xp <- list();  # deflated data sets
  cen <- list(); sc <- list();  # centering and scaling
  dx <- numeric(m)  # dimensionality of data domain
  for (mm in 1:m) {
    X[[mm]] <- scale(as.matrix(x[[mm]]), 
                     if (is.list(center)) center[[mm]] else center, 
                     if (is.list(scale_)) scale_[[mm]] else scale_
    )
    dx[mm] <- ncol(X[[mm]])
    
    cent <- attr(X[[mm]], "scaled:center")
    cen[[mm]] <- if(is.null(cent)) rep.int(0, dx[mm]) else cent  # follows cancor convention
    
    scal <- attr(X[[mm]], "scaled:scale")
    if(any(scal == 0))
      stop(paste("cannot rescale a constant column to unit variance in domain", mm))
    sc[[mm]] <- if(is.null(scal)) FALSE else scal
    
    Xp[[mm]] <- X[[mm]]
    attr(Xp[[mm]], "scaled:center") <- NULL
    attr(Xp[[mm]], "scaled:scale") <- NULL
  }
  
  corr <- array(NA, dim = c(m, m, nvar))  # additional explained correlation
  
  for (pp in seq(nvar)) {

    XpW <- matrix(NA, n, m)  
    
    for (mm in 1:m) {
      w <- W[[mm]][ , pp]
      XpW[ , mm] <- Xp[[mm]]%*%w
      
      # deflate data matrix
      q <- t(Xp[[mm]])%*%(X[[mm]]%*%w)
      q <- q/normv(q)
      Xp[[mm]] <- Xp[[mm]] - Xp[[mm]]%*%q%*%t(q)
    }
  
    corr[ , , pp] <- cor(XpW, XpW)
  }
    
  mcc <- list(cor = corr, coef = coef, center = cen, scale = sc, xp = Xp)
  class(mcc) <- "mcancor"
  return(mcc)
}