#  Copyright 2013, 2014 Christian Sigg
#  Copyright 1995-2012 The R Core Team
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

#' Non-Negative and Sparse Cumulative PCA
#' 
#' Performs a PCA-like analysis on the given data matrix, where
#' non-negativity and/or sparsity constraints are enforced on the principal axes
#' (PAs). In contrast to regular PCA, which greedily maximises the variance of
#' each principal component (PC), \code{nscumcomp} \emph{jointly} optimizes the
#' components such that the cumulative variance of all PCs is maximal.
#' 
#' \code{nscumcomp} computes all PCs jointly using expectation-maximization (EM)
#' iterations. The M-step is equivalent to minimizing the objective function
#' 
#' \deqn{\left\Vert \mathbf{X}-\mathbf{Z}\mathbf{W}^{\top}\right\Vert
#' _{F}^{2}+\gamma\left\Vert \mathbf{W}^{\top}\mathbf{W}-\mathbf{I}\right\Vert
#' _{F}^{2}}{norm(X - Z*t(W), "F")^2 + gamma*norm(t(W)*W - I, "F")^2}
#' 
#' w.r.t. the pseudo-rotation matrix \eqn{\mathbf{W}}{W}, where 
#' \eqn{\mathbf{Z}=\mathbf{X}\mathbf{W}\left(\mathbf{W}^\top\mathbf{W}\right)^{-1}}{Z=X*W*solve(t(W)*W)}
#' is the scores matrix modified to account for the non-orthogonality of 
#' \eqn{\mathbf{W}}{W}, \eqn{\mathbf{I}}{I} is the identity matrix and 
#' \code{gamma} is the Lagrange parameter associated with the ortho-normality 
#' penalty on \eqn{\mathbf{W}}{W}. Non-negativity of the loadings is achieved by
#' enforcing a zero lower bound in the L-BFGS-B algorithm used for the
#' minimization of the objective, and sparsity is achieved by a subsequent soft
#' thresholding of \eqn{\mathbf{W}}{W}.
#' 
#' @export nscumcomp
#' @param x a numeric matrix or data frame which provides the data for the
#'   analysis.
#' @param ... arguments passed to or from other methods.
nscumcomp <- function (x, ...) UseMethod("nscumcomp") 

#' @method nscumcomp default
#' @S3method nscumcomp default
#' @rdname nscumcomp
#' @param ncomp the number of principal components (PCs) to be computed. The 
#'   default is to compute a full basis for \code{x}.
#' @param omega a vector with as many entries as there are data samples, to 
#'   perform weighted PCA (analogous to weighted least-squares regression). The 
#'   default is an equal weighting of all samples.
#' @param k an upper bound on the total number of non-zero loadings of the 
#'   pseudo-rotation matrix \eqn{\mathbf{W}}{W}. \code{k} is increased if
#'   necessary to ensure at least one non-zero coefficient per principal axis.
#' @param nneg a logical value indicating whether the loadings should be 
#'   non-negative, i.e. the PAs should be constrained to the non-negative 
#'   orthant.
#' @param gamma a non-negative penalty on the divergence from orthonormality of 
#'   the pseudo-rotation matrix. The default is not to penalize, but a positive 
#'   value is sometimes necessary to avoid PAs collapsing onto each other.
#' @param center a logical value indicating whether the empirical mean of (the 
#'   columns of) \code{x} should be subtracted. Alternatively, a vector of 
#'   length equal to the number of columns of \code{x} can be supplied. The 
#'   value is passed to \code{\link{scale}}.
#' @param scale. a logical value indicating whether the columns of \code{x} 
#'   should be scaled to have unit variance before the analysis takes place. The
#'   default is \code{FALSE} for consistency with \code{prcomp}. Alternatively, 
#'   a vector of length equal to the number of columns of \code{x} can be 
#'   supplied. The value is passed to \code{\link{scale}}.
#' @param nrestart the number of random restarts for computing the 
#'   pseudo-rotation matrix via expectation-maximization (EM) iterations. The 
#'   solution achieving the minimum of the objective function over all random 
#'   restarts is kept. A value greater than one can help to avoid poor local 
#'   minima.
#' @param em_tol If the relative change of the objective is less than 
#'   \code{em_tol} between iterations, the EM procedure is asssumed to have 
#'   converged to a local optimum.
#' @param em_maxiter the maximum number of EM iterations to be performed. The EM
#'   procedure is terminated if either the \code{em_tol} or the 
#'   \code{em_maxiter} criterion is satisfied.
#' @param verbosity an integer specifying the verbosity level. Greater values 
#'   result in more output, the default is to be quiet.
#'   
#' @return \code{nscumcomp} returns a list with class \code{(nsprcomp, prcomp)} 
#'   containing the following elements: \item{sdev}{the additional standard 
#'   deviation explained by each component, see \code{\link{asdev}}.} 
#'   \item{rotation}{the matrix of non-negative and/or sparse loadings, 
#'   containing the principal axes as columns.} \item{x}{the scores matrix 
#'   \eqn{\mathbf{XW}}{X*W} containing the principal components as columns 
#'   (after centering and scaling if requested)} \item{center, scale}{the 
#'   centering and scaling used, or \code{FALSE}} \item{xp}{the deflated data 
#'   matrix corresponding to \code{x}} \item{q}{an orthonormal basis for the 
#'   principal subspace}
#'   
#'   The components are returned in order of decreasing variance for 
#'   convenience.
#'   
#' @note The PCA terminology is not consistent across the literature. Given a 
#'   zero mean data matrix \eqn{\mathbf{X}}{X} (with observations as rows) and a
#'   basis \eqn{\mathbf{W}}{W} of the principal subspace, we define the scores 
#'   matrix as \eqn{\mathbf{Z}=\mathbf{XW}}{Z=X*W} which contains the principal 
#'   components as its columns. The columns of the pseudo-rotation matrix 
#'   \eqn{\mathbf{W}}{W} are called the principal axes, and the elements of 
#'   \eqn{\mathbf{W}}{W} are called the loadings.
#'   
#' @seealso \code{\link{asdev}},  \code{\link{peav}}, \code{\link{nsprcomp}}, 
#'   \code{\link{scale}}
#'   
#' @example inst/nscumcomp_examples.R
nscumcomp.default <-
    function(x, ncomp = min(dim(x)), 
             omega = rep(1, nrow(x)),
             k = d*ncomp, nneg = FALSE, gamma = 0,
             center = TRUE, scale. = FALSE,
             nrestart = 5, em_tol = 1e-3, em_maxiter = 20,
             verbosity = 0, ...) 
{       
    d <- ncol(x); n <- nrow(x)
    if (k > d*ncomp) {
        k = d*ncomp
    }
    
    X <- as.matrix(x)
    X <- scale(X, center = center, scale = scale.)
    cen <- attr(X, "scaled:center")
    sc <- attr(X, "scaled:scale")
    if(any(sc == 0))
        stop("cannot rescale a constant/zero column to unit variance")
    
    obj_opt <- Inf
    for (rr in seq(nrestart)) {        
        res <- emcumca(X, omega, ncomp, k, nneg, gamma, NULL, em_tol, 
                       em_maxiter, verbosity)
        
        # variational renormalization
        S <- abs(res$W) > 0  # support of W
        if (any(!S))
            res <- emcumca(X, omega, ncomp, k, nneg, gamma, S, em_tol, 
                           em_maxiter, verbosity)
        
        # keep solution with minimum objective
        if (res$obj < obj_opt) {
            obj_opt <- res$obj
            W_opt <- res$W
        }
        if (verbosity > 0) {
            print(paste("minimum objective is ", format(obj_opt, digits = 4),
                        " at random restart ", rr-1, sep = ""))
        }
    }
    
    # sort components according to their additional variances
    sdev <- asdev(X, W_opt)$sdev
    if (length(sdev) > 1) {
        srt = sort(sdev, decreasing = TRUE, index.return = TRUE)
        sdev = srt$x
        W_opt <- W_opt[, srt$ix, drop=FALSE]
    }
    
    dimnames(W_opt) <-
        list(colnames(X), paste0("PC", seq_len(ncol(W_opt))))
    
    r <- list(sdev = sdev, rotation = W_opt,
              center = if(is.null(cen)) FALSE else cen,
              scale = if(is.null(sc)) FALSE else sc,
              x = X %*% W_opt)
    
    class(r) <- c("nsprcomp", "prcomp")
    return(r)
}

emcumca <- function(X, omega, ncomp, k, nneg, gamma, S, em_tol, em_maxiter, 
                    verbosity) {
    d <- ncol(X); n <- nrow(X)
    
    W_init <- function(W_) {  # if S is not null, W_ only contains the elements corresponding to the support
        if (is.null(S)) {
            W <- matrix(W_, d, ncomp)
        } else {
            W <- matrix(0, d, ncomp)   
            W[S] <- W_
        }
        return(W)
    }
    
    W <- matrix(rnorm(d*ncomp), d)
    if (nneg) {
        W <- abs(W)
    }
    if (!is.null(S)) {
        W[!S] <- 0
    }
    W <- W/t(matrix(rep(sqrt(colSums(W^2)), d), ncomp))
    
    obj_old <- Inf
    ii <- 0
    while(ii < em_maxiter) {   
        Z <- tryCatch(X%*%W%*%solve(t(W)%*%W), 
                 error = function(e) {
                     stop("Co-linear principal axes, try increasing the orthonormality penalty 'gamma'.")
                 }
        )  
        
        # objective function
        if (isTRUE(all.equal(omega, rep(1, nrow(X))))) {
            fn <- function(W_) {
                W <- W_init(W_)
                obj <- sum((X-Z%*%t(W))^2)
                if (gamma > 0)
                    obj <- obj + gamma*sum((t(W)%*%W - diag(ncomp))^2)
                return(obj)
            }   
        } else {
            fn <- function(W_) {
                W <- W_init(W_)
                obj <- omega%*%rowSums((X-Z%*%t(W))^2)
                if (gamma > 0)
                    obj <- obj + gamma*sum((t(W)%*%W - diag(ncomp))^2)
                return(obj)
            }  
        }
        
        # gradient
        if (isTRUE(all.equal(omega, rep(1, nrow(X))))) {
            tZZ <- t(Z)%*%Z
            tXZ <- t(X)%*%Z
            gr <- function(W_) {
                W <- W_init(W_)
                G <- 2*(W%*%tZZ - tXZ)
                if (gamma > 0)
                    G <- G + 4*gamma*(W%*%(t(W)%*%W) - W)
                if (is.null(S)) {
                    return(G)
                } else {
                    return(G[S])
                }
            }
        } else {
            OZ <- matrix(rep(omega, ncomp), n)*Z
            tZOZ <- t(Z)%*%OZ
            tXOZ <- t(X)%*%OZ
            gr <- function(W_) {
                W <- W_init(W_)
                G <- 2*(W%*%tZOZ - tXOZ)
                if (gamma > 0)
                    G <- G + 4*gamma*(W%*%(t(W)%*%W) - W)
                if (is.null(S)) {
                    return(G)
                } else {
                    return(G[S])
                }
            }
        }
        
        if (is.null(S)) {
            obj <- fn(W)     
        } else {
            obj <- fn(W[S])   
        }
        if ((obj == 0) || (abs(obj-obj_old)/obj < em_tol)) {
            break
        }
        obj_old <- obj
        
        if (verbosity > 1) {
            print(paste("approximation error: ", 
                        format(omega%*%rowSums((X-Z%*%t(W))^2), digits = 4),
                        " - ortho penalty: ", 
                        format(sum((t(W)%*%W - diag(ncomp))^2), digits = 4),
                        sep = "")
            )
        }
        
        control <- list()
        if (verbosity > 2) {
            control$trace <- 3    
        } else {
            control$trace <- 0
        }
        if (is.null(S)) {
            W_star <- optim(W, fn, gr, method = "L-BFGS-B",
                            lower = if (nneg) 0 else -Inf, control = control)$par    
        } else {
            W_star <- matrix(0, d, ncomp)
            W_star[S] <- optim(W[S], fn, gr, method = "L-BFGS-B",
                            lower = if (nneg) 0 else -Inf, control = control)$par
        }
        
        # soft thresholding
        if (is.null(S) && (k < d*ncomp)) {        
            aWstar <- abs(W_star)
            srt <- sort(aWstar, decreasing = TRUE)
            
            # increase k (by decreasing the threshold) if necessary to ensure 
            # at least one non-zero coefficient remains in each principal axis
            colmax <- apply(aWstar, 2, max)
            thr <- min(srt[k+1], 0.999*colmax)
            W <- pmax(aWstar-thr, 0)*sign(W_star)      
        } else {
            W <- W_star
        }
        
        norms <- sqrt(colSums(W^2))
        W <- W/t(matrix(rep(norms, d), ncomp))
        ii <- ii + 1
    }
    if ((verbosity > 0) && (ii == em_maxiter))
        print("maximum number of EM iterations reached before convergence")
    
    return(list(W=W, obj=obj))
}

#' @method nscumcomp formula
#' @S3method nscumcomp formula
#' @rdname nscumcomp
#' @param formula a formula with no response variable, referring only to numeric variables.
#' @param data an optional data frame (or similar: see
#'   \code{\link{model.frame}}) containing the variables in the
#'   formula \code{formula}.  By default the variables are taken from
#'   \code{environment(formula)}.
#' @param subset an optional vector used to select rows (observations) of the
#'   data matrix \code{x}.
#' @param na.action a function which indicates what should happen
#'   when the data contain \code{NA}s.  The default is set by
#'   the \code{na.action} setting of \code{\link{options}}, and is
#'   \code{\link{na.fail}} if that is unset. 
nscumcomp.formula <- function (formula, data = NULL, subset, na.action, ...) {
    mt <- terms(formula, data = data)
    if (attr(mt, "response") > 0L)
        stop("response not allowed in formula")
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    mf$... <- NULL
    mf[[1L]] <- as.name("model.frame")
    mf <- eval.parent(mf)
    ## this is not a `standard' model-fitting function,
    ## so no need to consider contrasts or levels
    if (.check_vars_numeric(mf))
        stop("PCA applies only to numerical variables")
    na.act <- attr(mf, "na.action")
    mt <- attr(mf, "terms")
    attr(mt, "intercept") <- 0L
    x <- model.matrix(mt, mf)
    res <- nscumcomp.default(x, ...)
    ## fix up call to refer to the generic, but leave arg name as `formula'
    cl[[1L]] <- as.name("nscumcomp")
    res$call <- cl
    if (!is.null(na.act)) {
        res$na.action <- na.act
        if (!is.null(sc <- res$x))
            res$x <- napredict(na.act, sc)
    }
    res
}
