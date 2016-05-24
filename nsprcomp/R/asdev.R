#  Copyright 2013, 2014 Christian Sigg
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

#' Additional Explained Standard Deviation
#' 
#' \code{asdev} computes the \emph{additional} standard deviation explained by each 
#' principal component, taking into account the possible non-orthogonality of
#' the pseudo-rotation matrix \eqn{\mathbf{W}}{W}.
#' 
#' The additional standard deviation of a component is measured after projecting
#' the corresponding principal axis to the ortho-complement space spanned by the
#' previous principal axes. This procedure ensures that the variance explained 
#' by non-orthogonal principal axes is not counted multiple times. If the 
#' principal axes are pairwise orthogonal (e.g. computed using standard PCA), 
#' the additional standard deviations are identical to the standard deviations 
#' of the columns of the scores matrix \eqn{\mathbf{XW}}{X*W}.
#' 
#' \code{asdev} is also useful to build a partial PCA model from
#' \eqn{\mathbf{W}}{W}, to be completed with additional components computed
#' using \code{\link{nsprcomp}}.
#' 
#' @export
#' @param x a numeric data matrix with the observations as rows
#' @param w a numeric data matrix with the principal axes as columns
#' @param center a logical value indicating whether the empirical mean of 
#'   \code{x} should be subtracted. Alternatively, a vector of length equal to 
#'   the number of columns of \code{x} can be supplied. The value is passed to 
#'   \code{\link{scale}}.
#' @param scale. a logical value indicating whether the columns of \code{x} 
#'   should be scaled to have unit variance before the analysis takes place. The
#'   default is \code{FALSE} for consistency with \code{prcomp}. Alternatively, 
#'   a vector of length equal to the number of columns of \code{x} can be 
#'   supplied.  The value is passed to \code{\link{scale}}.
#'   
#' @return \code{asdev} returns a list with class \code{(nsprcomp, prcomp)} 
#'   containing the following elements: \item{sdev}{the additional standard 
#'   deviation explained by each component} \item{rotation}{copied from the 
#'   input argument \code{w}} \item{x}{the scores matrix \eqn{\mathbf{XW}}{X*W}, 
#'   containing the principal components as columns (after centering and scaling
#'   if requested)} \item{center, scale.}{the centering and scaling used} 
#'   \item{xp}{the deflated data matrix corresponding to \code{x}} \item{q}{an 
#'   orthonormal basis for the principal subspace}
#'   
#' @note The PCA terminology is not consistent across the literature. Given a 
#'   zero mean data matrix \eqn{\mathbf{X}}{X} (with observations as rows) and a
#'   basis \eqn{\mathbf{W}}{W} of the principal subspace, we define the scores 
#'   matrix as \eqn{\mathbf{Z}=\mathbf{XW}}{Z=X*W} which contains the principal 
#'   components as its columns. The columns of the pseudo-rotation matrix 
#'   \eqn{\mathbf{W}}{W} are called the principal axes, and the elements of 
#'   \eqn{\mathbf{W}}{W} are called the loadings.
#'   
#' @references Mackey, L. (2009) Deflation Methods for Sparse PCA. In 
#'   \emph{Advances in Neural Information Processing Systems} (pp. 1017--1024).
asdev <- function(x, w, center = TRUE, scale. = FALSE) {
    
    X <- scale(x, center, scale.)
    cen <- attr(X, "scaled:center")
    sc <- attr(X, "scaled:scale")
    if(any(sc == 0))
        stop("cannot rescale a constant/zero column to unit variance")
    
    W <- w
    Q <- qr.Q(qr(W))  # orthonormal basis for the principal subspace
    Xp <- X  # deflated data matrix
    attr(Xp, "scaled:center") <- NULL
    attr(Xp, "scaled:scale") <- NULL

    nc <- ncol(W)
    sdev <- numeric(nc)
    for (cc in seq_len(nc)) {
        sdev[cc] <- sd(as.vector(Xp%*%W[ , cc]))  # explicit casting to avoid warning in old R versions
        Xp <- Xp - Xp%*%Q[ , cc]%*%t(Q[ , cc])   
    }
    
    nspc <- list(sdev = sdev, rotation = W,
                 center = if(is.null(cen)) FALSE else cen,
                 scale = if(is.null(sc)) FALSE else sc,
                 x = X%*%W, xp = Xp, q = Q)
    class(nspc) <- c("nsprcomp", "prcomp")
    return(nspc)
}