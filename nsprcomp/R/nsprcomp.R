#  Copyright 2012, 2013 Christian Sigg
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

#' Non-Negative and Sparse PCA
#' 
#' Performs a constrained principal component analysis,
#' where non-negativity and/or sparsity constraints are enforced on the principal axes.
#' The result is returned as an object of class \code{nsprcomp}, which inherits from
#' \code{prcomp}.
#' 
#' \code{nsprcomp} computes a principal component (PC) using expectation-maximization
#' iterations, where non-negativity of the loadings is achieved by projecting
#' the principal axis (PA) into the non-negative orthant, and sparsity of the 
#' loadings is achieved by soft thresholding (Sigg and Buhmann, 2008). 
#' 
#' Because constrained principal axes no longer correspond to true eigenvectors 
#' of the covariance matrix and are usually not pairwise orthogonal, special
#' attention needs to be paid when computing more than a single PC. The
#' algorithm implements the generalized deflation method proposed by
#' Mackey (2009) to maximize the additional variance of each
#' component. Given a basis of the space spanned by the previous PAs, the
#' variance of the PC is maximized after projecting the current PA to the
#' ortho-complement space of the basis. This procedure maximizes the
#' additional variance not explained by previous components, and is
#' identical to standard PCA if no sparsity or non-negativity constraints
#' are enforced on the PAs.
#' 
#' See the references for further details. 
#' 
#' @export nsprcomp
#' @param x a numeric matrix or data frame which provides the data 
#'   for the principal component analysis.
#' @param ... arguments passed to or from other methods.
nsprcomp <- function (x, ...) UseMethod("nsprcomp") 

#' @method nsprcomp default
#' @S3method nsprcomp default
#' @rdname nsprcomp
#' @param retx a logical value indicating whether the principal components, i.e.
#'   \code{x} projected into the principal subspace, should be returned.
#' @param ncomp the number of principal components (PCs) to be computed. With 
#'   the default setting, PCs are computed until \code{x} is fully deflated. 
#'   \code{ncomp} can be specified implicitly if \code{k} is given as a vector.
#' @param omega a vector with as many entries as there are data samples, to 
#'   perform weighted PCA (analogous to weighted least-squares regression). The 
#'   default is an equal weighting of all samples.
#' @param k either a scalar or a vector of length \code{ncomp}, specifying the 
#'   upper bounds on the cardinalities of the principal axes (PAs).
#' @param nneg a logical value indicating whether the loadings should be 
#'   non-negative, i.e. the PAs should be constrained to the non-negative 
#'   orthant.
#' @param center a logical value indicating whether the empirical mean of (the 
#'   columns) of \code{x} should be subtracted. Alternatively, a vector of 
#'   length equal to the number of columns of \code{x} can be supplied. The 
#'   value is passed to \code{\link{scale}}.
#' @param scale. a logical value indicating whether the columns of \code{x} 
#'   should be scaled to have unit variance before the analysis takes place. The
#'   default is \code{FALSE} for consistency with \code{prcomp}. Alternatively, 
#'   a vector of length equal to the number of columns of \code{x} can be 
#'   supplied. The value is passed to \code{\link{scale}}.
#' @param tol a threshold indicating the magnitude below which components should
#'   be omitted. Components are omitted if their standard deviations are less 
#'   than or equal to \code{tol} times the standard deviation of the first 
#'   component. With the default \code{NULL} setting, no components are omitted.
#'   With \code{tol = 0} or \code{tol = sqrt(.Machine$double.eps)}, essentially 
#'   constant components are omitted.
#' @param nrestart the number of random restarts for computing the principal 
#'   component via expectation-maximization (EM) iterations. The solution 
#'   achieving maximum standard deviation over all random restarts is kept. A 
#'   value greater than one can help to avoid poor local maxima.
#' @param em_tol If the relative change of the objective is less than 
#'   \code{em_tol} between iterations, the EM procedure is asssumed to have 
#'   converged to a local optimum.
#' @param em_maxiter the maximum number of EM iterations to be performed. The EM
#'   procedure is terminated if either the \code{em_tol} or the 
#'   \code{em_maxiter} criterion is satisfied.
#' @param partial_model \code{NULL} or an object of class \code{nsprcomp}. The 
#'   computation can be continued from a partial model by providing a 
#'   \code{nsprcomp} object (either from a previous run of this function or from
#'   \code{\link{asdev}}) and setting \code{ncomp} to a value greater than the 
#'   number of components contained in the partial model. See the examples for 
#'   an illustration.
#' @param verbosity an integer specifying the verbosity level. Greater values 
#'   result in more output, the default is to be quiet.
#'   
#' @return \code{nsprcomp} returns a list with class \code{(nsprcomp, prcomp)} 
#'   containing the following elements: \item{sdev}{the additional standard 
#'   deviation explained by each component, see \code{\link{asdev}}.} 
#'   \item{rotation}{the matrix of non-negative and/or sparse loadings, 
#'   containing the principal axes as columns.} \item{x}{the scores matrix 
#'   \eqn{\mathbf{XW}}{X*W} containing the principal components as columns 
#'   (after centering and scaling if requested). For the formula method, 
#'   \code{\link{napredict}} is applied to handle the treatment of values 
#'   omitted by the \code{na.action}.} \item{center, scale}{the centering and 
#'   scaling used, or \code{FALSE}} \item{xp}{the deflated data matrix 
#'   corresponding to \code{x}} \item{q}{an orthonormal basis for the principal 
#'   subspace}
#'   
#' @note The PCA terminology is not consistent across the literature. Given a 
#'   zero mean data matrix \eqn{\mathbf{X}}{X} (with observations as rows) and a
#'   basis \eqn{\mathbf{W}}{W} of the principal subspace, we define the scores 
#'   matrix as \eqn{\mathbf{Z}=\mathbf{XW}}{Z=X*W} which contains the principal 
#'   components as its columns. The columns of the pseudo-rotation matrix
#'   \eqn{\mathbf{W}}{W} are called the principal axes, and the elements of
#'   \eqn{\mathbf{W}}{W} are called the loadings.
#'   
#'   Deflating the data matrix accumulates numerical errors over successive PCs.
#'   
#' @references Sigg, C. D. and Buhmann, J. M. (2008) Expectation-Maximization 
#'   for Sparse and Non-Negative PCA. In \emph{Proceedings of the 25th 
#'   International Conference on Machine Learning} (pp. 960--967).
#' @references Mackey, L. (2009) Deflation Methods for Sparse PCA. In 
#'   \emph{Advances in Neural Information Processing Systems} (pp. 1017--1024).
#'   
#' @seealso  \code{\link{asdev}},  \code{\link{peav}}, \code{\link{prcomp}}, 
#'   \code{\link{scale}}
#'   
#' @example inst/nsprcomp_examples.R
nsprcomp.default <- function(x, retx = TRUE, ncomp = min(dim(x)), 
                             omega = rep(1, nrow(x)),
                             k = ncol(x), nneg = FALSE, 
                             center = TRUE, scale. = FALSE, tol = NULL,
                             nrestart = 5, em_tol = 1e-3, em_maxiter = 100,
                             partial_model = NULL,
                             verbosity = 0, ...) {  
    
    d <- ncol(x); n <- nrow(x)
    if (length(k) == 1) {
        k <- rep(k, ncomp)
    } else {
        ncomp <- length(k)
    }
    
    X <- as.matrix(x)
    X <- scale(X, center = center, scale = scale.)
    cen <- attr(X, "scaled:center")
    sc <- attr(X, "scaled:scale")
    if(any(sc == 0))
        stop("cannot rescale a constant/zero column to unit variance")
    
    sdev <- rep(0, ncomp)  # additional explained standard deviations
    W <- matrix(0, d, ncomp)  # principal axes as columns
    Q <- matrix(0, d, ncomp)  # orthonormal basis for the principal subspace
    Xp <- X  # data matrix, to be deflated if ncomp > 1
    attr(Xp, "scaled:center") <- NULL
    attr(Xp, "scaled:scale") <- NULL
    ccs <- seq(ncomp)

    # start from partially completed model
    if (!is.null(partial_model)) {
        if (!("nsprcomp" %in% class(partial_model))) {
            stop("argument 'partial_model' must have class 'nsprcomp'")
        }

        pcomp <- length(partial_model$sdev)  # number of components in partial model
        if (pcomp >= ncomp)
            stop("'ncomp' must exceed the number of canonical variables in the partial model")
        sdev[1:pcomp] <- partial_model$sdev
        W[ , 1:pcomp] <- partial_model$rotation
        Q[ , 1:pcomp] <- partial_model$q
        Xp <- partial_model$xp
        ccs <- seq(pcomp+1, ncomp)
    }

    for (cc in ccs) {
        obj_opt <- -Inf
        for (rr in seq(nrestart)) {
            res <- empca(Xp, Q[ , seq_len(cc-1), drop=FALSE], omega, k[cc], 
                         nneg, em_tol, em_maxiter, verbosity)
            
            # variational renormalization
            supp <- abs(res$w) > 0  # support of w
            w <- rep(0, d)
            res <- empca(Xp[ , supp, drop=FALSE], 
                         Q[supp, seq_len(cc-1), drop=FALSE], 
                         omega, sum(supp), nneg, em_tol, em_maxiter, 
                         verbosity)
            w[supp] <- res$w                
            
            # keep solution with maximum objective
            obj_cur <- res$obj
            if (obj_cur > obj_opt) {
                obj_opt <- obj_cur
                w_opt <- w
            }
            if (verbosity > 0) {
                print(paste("component ", cc, ": ",
                            "maximum objective is ", format(obj_opt, digits = 4),
                            " at random restart ", rr-1, sep = ""))
            }
        }
        w <- w_opt  
        W[ ,cc] <- w
        sdev[cc] <- sd(as.vector(Xp%*%w))  # explicit casting to avoid warning in old R versions
        
        # early stopping at previous component due to standard deviation threshold
        if (!is.null(tol) && sdev[cc] < tol*sdev[1]) {
            if (verbosity > 0) {
                print("component magnitude below threshold, less than 'ncomp' components could be computed")
            }
            sdev <- sdev[1:(cc-1)]
            W <- W[ , 1:(cc-1), drop=FALSE]
            Q <- Q[ , 1:(cc-1), drop=FALSE]
            break
        }

        # deflate data matrix
        if (cc > 1) {
            q <- w - Q[ , 1:(cc-1)]%*%(t(Q[ , 1:(cc-1)])%*%w) 
        } else {
            q <- w
        }
        q <- q/normv(q)
        Q[ , cc] <- q
        Xp <- Xp - Xp%*%q%*%t(q)
        
        # early stopping at current component due to fully deflated data matrix
        if (cc < ncomp && all(abs(Xp) < 1e-14)) {
            if (verbosity > 0) {
                print("data matrix is fully deflated, less than 'ncomp' components could be computed")            
            }
            sdev <- sdev[1:cc]
            W <- W[ , 1:cc, drop=FALSE]
            Q <- Q[ , 1:cc, drop=FALSE]
            break
        }
    }
    
    dimnames(W) <- list(colnames(X), paste0("PC", seq_len(ncol(W))))
    nspc <- list(sdev = sdev, rotation = W,
              center = if(is.null(cen)) FALSE else cen,
              scale = if(is.null(sc)) FALSE else sc,
              xp = Xp, q = Q)
    if (retx) nspc$x <- X%*%W
    class(nspc) <- c("nsprcomp", "prcomp")
    return(nspc)
}

empca <- function(Xp, Q, omega, k, nneg, em_tol, em_maxiter, verbosity) {
    
    d <- ncol(Xp); n <- nrow(Xp)
    
    w <- rnorm(d); 
    if (nneg) {
        w <- abs(w)
    }
    w <- w/normv(w)
    
    obj_old <- -Inf
    ii <- 0
    while(ii < em_maxiter) {  
        z <- Xp%*%w
        
        obj <- sum(z^2)
        if (abs(obj - obj_old)/obj < em_tol) {
            break
        }
        obj_old <- obj
        
        w_star <- as.vector(t(Xp)%*%(omega*z)/as.vector(t(z)%*%(omega*z)))
        if (nneg) {
            w_star[w_star < 0] <- 0
        }
        if (k<d) {
            # soft thresholding
            s <- sort(abs(w_star), decreasing = TRUE)
            w <- pmax(abs(w_star) - s[k+1], 0)*sign(w_star)  
        } else {
            w <- w_star
        }
        
        if (ncol(Q) > 0) {
            wo <- t(Q) %*% w
            w <- w/sqrt(t(w)%*%w - t(wo)%*%wo)    
        } else {
            w <- w/normv(w)
        }
        
        ii <- ii + 1
    }
    if ((verbosity > 0) && (ii == em_maxiter))
        print("maximum number of EM iterations reached before convergence")
    
    w <- w/normv(w)
    return(list(w=w, obj=obj))
}

#' @method nsprcomp formula
#' @S3method nsprcomp formula
#' @rdname nsprcomp
#' @param formula a formula with no response variable, referring only to numeric variables.
#' @param data an optional data frame (or similar: see
#'   \code{\link{model.frame}}) containing the variables in the
#'   \code{formula}. By default the variables are taken from
#'   \code{environment(formula)}.
#' @param subset an optional vector used to select rows (observations) of the
#'   data matrix \code{x}.
#' @param na.action a function which indicates what should happen
#'   when the data contain \code{NA}s.  The default is set by
#'   the \code{na.action} setting of \code{\link{options}}, and is
#'   \code{\link{na.fail}} if that is unset. The \sQuote{factory-fresh}
#'   default is \code{\link{na.omit}}.
nsprcomp.formula <- function (formula, data = NULL, subset, na.action, ...) {
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
    res <- nsprcomp.default(x, ...)
    ## fix up call to refer to the generic, but leave arg name as `formula'
    cl[[1L]] <- as.name("nsprcomp")
    res$call <- cl
    if (!is.null(na.act)) {
        res$na.action <- na.act
        if (!is.null(sc <- res$x))
            res$x <- napredict(na.act, sc)
    }
    res
}

# copied as is from R 3.0.1 src/library/stats/R/prcomp.R to avoid a warning
# when importing using :::
.check_vars_numeric <- function(mf)
{
    ## we need to test just the columns which are actually used.
    mt <- attr(mf, "terms")
    mterms <- attr(mt, "factors")
    mterms <- rownames(mterms)[apply(mterms, 1L, function(x) any(x > 0L))]
    any(sapply(mterms, function(x) is.factor(mf[,x]) || !is.numeric(mf[,x])))
}