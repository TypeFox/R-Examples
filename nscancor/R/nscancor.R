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

#' Non-Negative and Sparse CCA
#' 
#' Performs a canonical correlation analysis (CCA) where constraints such as 
#' non-negativity or  sparsity are enforced on the canonical vectors. The result
#' of the analysis is returned as a list of class \code{nscancor}, which 
#' contains a superset of the elements returned by \code{\link{cancor}}.
#' 
#' \code{nscancor} computes the canonical vectors (called \code{xcoef} and 
#' \code{ycoef}) using iterated regression steps, where the constraints suitable
#' for each domain are enforced by choosing the avvropriate regression method. 
#' See Sigg et al. (2007) for an early avvlication of the principle (not yet 
#' including generalized deflation).
#' 
#' Because constrained canonical vectors no longer correspond to true 
#' eigenvectors of the cross-covariance matrix and are usually not pairwise 
#' conjugate (i.e. the canonical variables are not uncorrelated), special 
#' attention needs to be paid when computing more than a single pair of 
#' canonical vectors. \code{nscancor} implements a generalized deflation (GD) 
#' scheme which builds on GD for PCA as proposed by Mackey (2009). For each 
#' domain, a basis of the space spanned by the previous canonical variables is 
#' computed. Then, the correlation of the current pair of canonical variables is
#' maximized after projecting each current canonical vector to the 
#' ortho-complement space of its respective basis. This procedure maximizes the 
#' additional correlation not explained by previous canonical variables, and is 
#' identical to standard CCA if the canonical vectors are the eigenvectors of 
#' the cross-covariance matrix.
#' 
#' See the references for further details.
#' 
#' @export nscancor
#' @param x a numeric matrix which provides the data from the first domain
#' @param y a numeric matrix which provides the data from the second domain
#' @param xcenter a logical value indicating whether the empirical mean of (each
#'   column of) \code{x} should be subtracted. Alternatively, a vector of length
#'   equal to the number of columns of \code{x} can be supplied. The value is 
#'   passed to \code{\link{scale}}.
#' @param ycenter analogous to \code{xcenter}
#' @param xscale a logical value indicating whether the columns of \code{x} 
#'   should be scaled to have unit variance before the analysis takes place. The
#'   default is \code{FALSE} for consistency with \code{cancor}. Alternatively, 
#'   a vector of length equal to the number of columns of \code{x} can be 
#'   supplied. The value is passed to \code{\link{scale}}.
#' @param yscale analogous to \code{xscale}
#' @param nvar the number of canonical variables to be computed for each domain.
#'   With the default setting, canonical variables are computed until either 
#'   \code{x} or \code{y} is fully deflated.
#' @param xpredict the regression function to predict the canonical variable for
#'   \code{x}, given \code{y}. The formal arguments are the design matrix 
#'   \code{y}, the regression target \code{xc} as the current canonical variable
#'   for \code{x}, and \code{cc} as a counter of the current pair of canonical 
#'   variables (e.g. for enforcing different constraints for different canonical
#'   vectors). See the examples for an illustration.
#' @param ypredict analogous to \code{xpredict}
#' @param cor_tol a threshold indicating the magnitude below which canonical 
#'   variables should be omitted. Variables are omitted if their explained 
#'   correlations are less than or equal to \code{cor_tol} times the correlation
#'   of the first pair of canonical variables. With the default \code{NULL} 
#'   setting, no variables are omitted.
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
#' @param partial_model \code{NULL} or an object of class \code{nscancor}. The 
#'   computation can be continued from a partial model by providing an 
#'   \code{nscancor} object (either from a previous run of this function or from
#'   \code{\link{acor}}) and setting \code{nvar} to a value greater than the 
#'   number of canonical variables contained in the partial model. See the
#'   examples for an illustration.
#' @param verbosity an integer specifying the verbosity level. Greater values 
#'   result in more output, the default is to be quiet.
#'   
#' @return \code{nscancor} returns a list of class \code{nscancor} containing 
#'   the following elements: \item{cor}{the additional correlation explained by 
#'   each pair of canonical variables, see \code{\link{acor}}.} \item{xcoef}{the
#'   matrix containing the canonical vectors related to \code{x} as its columns}
#'   \item{ycoef}{analogous to \code{xcoef}} \item{xcenter}{if \code{xcenter} is
#'   \code{TRUE} the centering vector, else the zero vector (in accordance with 
#'   \code{cancor})} \item{ycenter}{analogous to \code{xcenter}} 
#'   \item{xscale}{if \code{xscale} is \code{TRUE} the scaling vector, else 
#'   FALSE } \item{yscale}{analogous to \code{xscale}} \item{xp}{the deflated 
#'   data matrix corresponding to \code{x}} \item{yp}{anologous to \code{xp}}
#'   
#' @references Sigg, C. and Fischer, B. and Ommer, B. and Roth, V. and Buhmann, 
#'   J. (2007) Nonnegative CCA for Audiovisual Source Separation. In 
#'   \emph{Proceedings of the 2007 IEEE Workshop on Machine Learning for Signal 
#'   Processing} (vv. 253--258).
#' @references Mackey, L. (2009) Deflation Methods for Sparse PCA. In 
#'   \emph{Advances in Neural Information Processing Systems} (vv. 1017--1024).
#'   
#' @seealso  \code{\link{acor}}, \code{\link{cancor}}, \code{\link{scale}}
#'   
#' @example inst/nscancor_examples.R
nscancor <- function(x, y, xcenter = TRUE, ycenter = TRUE, 
                     xscale = FALSE, yscale = FALSE, nvar = min(dim(x), dim(y)),
                     xpredict, ypredict, 
                     cor_tol = NULL, nrestart = 10, iter_tol = 1e-3, iter_max = 30,
                     partial_model = NULL,
                     verbosity = 0) {  
  
    n <- nrow(x)
    dx <- ncol(x)
    dy <- ncol(y)
    
    X <- as.matrix(x)
    X <- scale(X, center = xcenter, scale = xscale)
    xcen <- attr(X, "scaled:center")
    xsc <- attr(X, "scaled:scale")
    Y <- as.matrix(y)
    Y <- scale(Y, center = ycenter, scale = yscale)
    ycen <- attr(Y, "scaled:center")
    ysc <- attr(Y, "scaled:scale")
    if(any(xsc == 0) || any(ysc == 0))
        stop("cannot rescale a constant/zero column to unit variance")
    
    corr <- rep(0, nvar)  # additional explained correlation
    W <- matrix(0, dx, nvar)  # canonical vectors for X
    V <- matrix(0, dy, nvar)  # canonical vectors for Y
    Xp <- X  # deflated X
    attr(Xp, "scaled:center") <- NULL
    attr(Xp, "scaled:scale") <- NULL
    Yp <- Y  # deflated Y
    attr(Yp, "scaled:center") <- NULL
    attr(Yp, "scaled:scale") <- NULL
    ccs <- seq(nvar)

    # start from partially completed model
    if (!is.null(partial_model)) {
      if (class(partial_model) != "nscancor") {
        stop("argument 'partial_model' must be of class 'nscancor'")
      }
      
      pvar <- length(partial_model$cor)  # number of canonical variables in partial model
      if (pvar >= nvar)
        stop("'nvar' must exceed the number of canonical variables in the partial model")
      corr[1:pvar] <- partial_model$cor
      W[ , 1:pvar] <- partial_model$xcoef
      V[ , 1:pvar] <- partial_model$ycoef
      Xp <- partial_model$xp
      Yp <- partial_model$yp
      ccs <- seq(pvar+1, nvar)
    } 
    
    for (cc in ccs) {
        obj_opt <- -Inf
        for (rr in seq(nrestart)) {
            
            res <- nscc_inner(X, Xp, Y, Yp, xpredict, ypredict, cc, iter_tol, 
                         iter_max, verbosity)               
            
            # keep solution with maximum objective
            if (res$obj > obj_opt) {
                obj_opt <- res$obj
                w_opt <- res$w
                v_opt <- res$v
            }
            if (verbosity > 0) {
                print(paste("canonical variable ", cc, ": ",
                            "maximum objective is ", format(obj_opt, digits = 4),
                            " at random restart ", rr-1, sep = ""))
            }
        }
        w <- w_opt  
        v <- v_opt
        W[ ,cc] <- w
        V[ ,cc] <- v
        corr[cc] <- cor(Xp%*%w, Yp%*%v)
        
        # early stopping at previous component due to correlation threshold
        if (!is.null(cor_tol) && corr[cc] < cor_tol*corr[1]) {
            corr <- corr[1:(cc-1)]
            W <- W[ , 1:(cc-1), drop=FALSE]
            V <- V[ , 1:(cc-1), drop=FALSE]
            break
        }    
        
        # deflate data matrices
        qx <- t(Xp)%*%(X%*%w)
        qx <- qx/normv(qx)
        Xp <- Xp - Xp%*%qx%*%t(qx)  
        
        qy <- t(Yp)%*%(Y%*%v)
        qy <- qy/normv(qy)
        Yp <- Yp - Yp%*%qy%*%t(qy)
        
        # early stopping at current component due to fully deflated data matrix
        if (cc < nvar && (all(abs(Xp) < 1e-14) || all(abs(Yp) < 1e-14))) { 
            if (verbosity > 0) {
                print("at least one data matrix is fully deflated, less than 'nvar' pairs of variables could be computed")            
            }
            corr <- corr[1:cc]
            W <- W[ , 1:cc, drop=FALSE]
            V <- V[ , 1:cc, drop=FALSE]
            break
        }
    }
    
    rownames(W) <- colnames(X)
    rownames(V) <- colnames(Y) 
    nscc <- list(cor = corr, xcoef = W, ycoef = V,
                 xcenter = if(is.null(xcen)) rep.int(0, dx) else xcen,  # return value follows cancor interface
                 xscale = if(is.null(xsc)) FALSE else xsc,
                 ycenter = if(is.null(ycen)) rep.int(0, dy) else ycen,
                 yscale = if(is.null(ysc)) FALSE else ysc,
                 xp = Xp, yp = Yp)
    class(nscc) <- "nscancor"
    return(nscc)
}

nscc_inner <- function(X, Xp, Y, Yp, xpredict, ypredict, cc, iter_tol, iter_max, 
                  verbosity) {
    
    n <- nrow(X)
    dx <- ncol(X)
    dy <- ncol(Y)
    
    # initialize w and v as random unit vectors, non-negative if the prediction
    # functions return non-negative vectors
    w <- rnorm(dx);
    v <- rnorm(dy);
    if (all(ypredict(Xp, Yp[ , 1], cc) >= 0)) {
        w <- abs(w)
    }
    if (all(xpredict(Yp, Xp[ , 1], cc) >= 0)) {
        v <- abs(v)
    }
    w <- w/normv(Xp%*%w)     
    v <- v/normv(Yp%*%v)     
    
    obj_old <- -Inf
    ii <- 0
    while(ii < iter_max) { 
        obj <- t(Xp%*%w)%*%(Yp%*%v)
        if (abs(obj - obj_old)/obj < iter_tol) {
            break
        }
        obj_old <- obj
        
        w <- ypredict(Xp, Yp%*%v, cc)
        if (all(w == 0))
            stop("w collapsed to the zero vector, try relaxing the constraints")
        w <- w/normv(Xp%*%w)
        
        v <- xpredict(Yp, Xp%*%w, cc)
        if (all(v == 0))
            stop("v collapsed to the zero vector, try relaxing the constraints")
        v <- v/normv(Yp%*%v)
        
        ii <- ii + 1
    }
    if ((verbosity > 0 ) && (ii == iter_max))
        print("maximum number of iterations reached before convergence")
    
    w <- w/normv(X%*%w)     
    v <- v/normv(Y%*%v) 
    return(list(obj=obj, w=w, v=v))
}