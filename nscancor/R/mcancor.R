#  Copyright 2013 Christian Sigg
#  Copyright 1995-2013 The R Core Team
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

#' Non-Negative and Sparse Multi-Domain CCA
#' 
#' Performs a canonical correlation analysis (CCA) on multiple data domains, 
#' where constraints such as non-negativity or sparsity are enforced on the 
#' canonical vectors. The result
#' of the analysis is returned as a list of class \code{mcancor}.
#' 
#' \code{mcancor} generalizes \code{\link{nscancor}} to the case where more than
#' two data domains are available for an analysis. Its objective is to maximize 
#' the sum of all pairwise correlations of the canonical variables.
#' 
#' @export
#' @param x a list of numeric matrices which contain the data from the different
#'   domains
#' @param center a list of logical values indicating whether the empirical mean 
#'   of (each column of) the corresponding data matrix should be subtracted. 
#'   Alternatively, a list of vectors can be supplied, where each vector 
#'   specifies the mean to be subtracted from the corresponding data matrix. 
#'   Each list element is passed to \code{\link{scale}}.
#' @param scale_ a list of logical values indicating whether the columns of the 
#'   corresponding data matrix should be scaled to have unit variance before the
#'   analysis takes place. The default is \code{FALSE} for consistency with 
#'   \code{nscancor}. Alternatively, a list of vectors can be supplied, where 
#'   each vector specifies the standard deviations used to rescale the columns 
#'   of the corresponding data matrix. Each list element is passed to 
#'   \code{\link{scale}}.
#' @param nvar the number of canonical variables to be computed for each domain.
#'   With the default setting, canonical variables are computed until at least 
#'   one data matrix is fully deflated.
#' @param predict a list of regression functions to predict the sum of the 
#'   canonical variables of all other domains. The formal arguments for each 
#'   regression function are the design matrix \code{x} corresponding to the 
#'   data from the current domain, the regression target \code{sc} as the sum of
#'   the canonical variables for all other domains, and \code{cc} as a counter 
#'   of which canonical variable is currently computed (e.g. for enforcing 
#'   different constraints for subsequent canonical vectors of a given domain). 
#'   See the examples for an illustration.
#' @param cor_tol a threshold indicating the magnitude below which canonical 
#'   variables should be omitted. Variables are omitted if the sum of all their 
#'   correlations are less than or equal to \code{cor_tol} times the sum of all 
#'   correlations of the first canonical variables of all domains. With the 
#'   default \code{NULL} setting, no variables are omitted.
#' @param nrestart the number of random restarts for computing the canonical 
#'   variables via iterated regression steps. The solution achieving maximum 
#'   explained correlation over all random restarts is kept. A value greater 
#'   than one can help to avoid poor local maxima.
#' @param iter_tol If the relative change of the objective is less than 
#'   \code{iter_tol} between iterations, the procedure is asssumed to have 
#'   converged to a local optimum.
#' @param iter_max the maximum number of iterations to be performed. The 
#'   procedure is terminated if either the \code{iter_tol} or the 
#'   \code{iter_max} criterion is satisfied.
#' @param partial_model \code{NULL} or an object of class \code{mcancor}. The 
#'   computation can be continued from a partial model by providing an 
#'   \code{mcancor} object (either from a previous run of this function or from
#'   \code{\link{macor}}) and setting \code{nvar} to a value greater than the 
#'   number of canonical variables contained in the partial model. See the
#'   examples for an illustration.
#' @param verbosity an integer specifying the verbosity level. Greater values 
#'   result in more output, the default is to be quiet.
#'   
#' @return \code{mcancor} returns a list of class \code{mcancor} with the following elements: 
#'   \item{cor}{a multi-dimensional array containing the additional correlations
#'   explained by each pair of canonical variables. The first two dimensions
#'   correspond to the domains, and the third dimension corresponds to the
#'   different canonical variables per domain (see also \code{\link{macor}}).} 
#'   \item{coef}{a list of matrices containing the canonical vectors related to 
#'   each data domain. The canonical vectors are stored as the columns of each 
#'   matrix.} \item{center}{the list of empirical means used to center the data
#'   matrices} \item{scale}{the list of empirical standard deviations used to
#'   scale the data matrices} \item{xp}{the list of deflated 
#'   data matrices corresponding to \code{x}}
#'   
#' @seealso  \code{\link{macor}}, \code{\link{nscancor}}, \code{\link{scale}}
#'   
#' @example inst/mcancor_examples.R
mcancor <- function (x, center = TRUE, scale_ = FALSE, nvar = min(sapply(x, dim)), 
                     predict, cor_tol = NULL, nrestart = 10, iter_tol = 1e-3, 
                     iter_max = 30, partial_model = NULL, verbosity = 0) {
  
  m <- length(x)  # number of domains
  
  X <- list();  # data sets
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
  
  n <- nrow(X[[1]])  # number of observations
  if (is.null(nvar))
    nvar <- min(n, dx)
  
  corr <- array(NA, dim = c(m, m, nvar))  # additional explained correlation
  W <- list()  # canonical vectors 
  for (mm in 1:m) {
    W[[mm]] <- matrix(NA, dx[mm], nvar)
    rownames(W[[mm]]) <- colnames(X[[mm]])
  }
  ccs <- seq(nvar)
  
  # start from partially completed model
  if (!is.null(partial_model)) {
    if (class(partial_model) != "mcancor") {
      stop("argument 'partial_model' must be of class 'mcancor'")
    }
    
    pvar <- dim(partial_model$cor)[3]  # number of canonical variables in partial model
    if (pvar >= nvar)
      stop("'nvar' must exceed the number of canonical variables in the partial model")
    corr[ , , 1:pvar] <- partial_model$cor
    
    for (mm in 1:m) {
      W[[mm]][ , 1:pvar] <- partial_model$coef[[mm]]
    }
    
    Xp <- partial_model$xp
    ccs <- seq(pvar+1, nvar)
  }
  
  for (cc in ccs) {
    obj_opt <- -Inf
    for (rr in seq(nrestart)) {
      
      res <- mcc_inner(X, Xp, dx, predict, cc, iter_tol, iter_max, verbosity)               
      
      # keep solution with maximum objective
      if (res$obj > obj_opt) {
        obj_opt <- res$obj
        w_opt <- res$w
        XpW <- res$XpW
      }
      if (verbosity > 0) {
        print(paste("canonical variable ", cc, ": ",
                    "maximum objective is ", format(obj_opt, digits = 4),
                    " at random restart ", rr-1, sep = ""))
      }
    }
    
    corr[ , , cc] <- cor(XpW, XpW)  
    
    # early stopping at previous component due to correlation threshold
    sum_corr_pp <- sum(corr[ , , cc] - diag(m))/2
    sum_corr_1 <- sum(corr[ , , 1] - diag(m))/2
    if (!is.null(cor_tol) && sum_corr_pp < cor_tol*sum_corr_1) {
        corr <- corr[ , , 1:(cc-1), drop=FALSE]
        for (mm in 1:m) {
            W[[mm]] <- W[[mm]][ , 1:(cc-1), drop=FALSE]  
        }
        break
    }    
    
    for (mm in 1:m) {
      w <- w_opt[[mm]]
      W[[mm]][ , cc] <- w
      
      # deflate data matrix
      q <- t(Xp[[mm]])%*%(X[[mm]]%*%w)
      q <- q/normv(q)
      Xp[[mm]] <- Xp[[mm]] - Xp[[mm]]%*%q%*%t(q) 
    }
    
    # early stopping at current component due to fully deflated data matrix
    deflated <- logical(m)
    for (mm in 1:m) {
      deflated[mm] <- all(abs(Xp[[mm]]) < 1e-14)
    }
    if (cc < nvar && any(deflated)) { 
      if (verbosity > 0) {
        print("at least one data matrix is fully deflated, less than 'nvar' pairs of variables could be computed")            
      }
      corr <- corr[ , , 1:cc, drop=FALSE]
      for (mm in 1:m) {
        W[[mm]] <- W[[mm]][ , 1:cc, drop=FALSE]  
      }
      break
    }
  }
  
  mcc <- list(cor = corr, coef = W, center = cen, scale = sc, xp = Xp)
  class(mcc) <- "mcancor"
  return(mcc)
}

mcc_inner <- function(X, Xp, dx, predict, cc, iter_tol, iter_max, verbosity) {
  
  m <- length(X)  # number of domains
  n <- nrow(X[[1]])  # number of observations
  
  # initialize w as a random unit vector, non-negative if the prediction
  # function returns a non-negative vector
  w <- list()
  for (mm in 1:m) {
    v <- rnorm(dx[mm]);  
    if (all(predict[[mm]](Xp[[mm]], Xp[[mm]][ , 1], cc) >= 0)) {
      v <- abs(v)
    }
    w[[mm]] <- v/normv(Xp[[mm]]%*%v)   
  }    
  
  obj_old <- -Inf
  ii <- 0
  XpW <- matrix(NA, n, m)
  while(ii < iter_max) { 
    
    for (mm in 1:m) {
      XpW[ , mm] <- Xp[[mm]]%*%w[[mm]]
    }
    obj <- sum(t(XpW)%*%(XpW))
    if (abs(obj - obj_old)/obj < iter_tol) {
      break
    }
    obj_old <- obj
    
    for (mm in 1:m) {
      v <- predict[[mm]](Xp[[mm]], rowSums(XpW[ , -mm, drop=FALSE]), cc)  
      if (all(v == 0))
        stop("w collapsed to the zero vector, try relaxing the constraints")
      
      w[[mm]] <- v/normv(Xp[[mm]]%*%v)
      XpW[ , mm] <- Xp[[mm]]%*%w[[mm]]
    }
    
    ii <- ii + 1
  }
  if ((verbosity > 0 ) && (ii == iter_max))
    print("maximum number of iterations reached before convergence")
  
  for (mm in 1:m) {
    w[[mm]] <- w[[mm]]/normv(X[[mm]]%*%w[[mm]])       
  }
  
  return(list(obj = obj, w = w, XpW = XpW))
}