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

#' Additional Explained Correlation
#' 
#' \code{acor} computes the additional standard correlation explained by each 
#' canonical variable, taking into account the possible non-conjugacy of the 
#' canonical vectors. The result of the analysis is returned as a list of class 
#' \code{nscancor}.
#' 
#' The additional correlation is measured after projecting the corresponding 
#' canonical vectors to the ortho-complement space spanned by the previous 
#' canonical variables. This procedure ensures that the correlation explained by
#' non-conjugate canonical vectors is not counted multiple times. See Mackey
#' (2009) for a presentation of generalized deflation in the context of
#' principal component analysis (PCA), which was adapted here to CCA.
#' 
#' \code{acor} is also useful to build a partial CCA model, to be completed with
#' additional canonical variables computed using \code{\link{nscancor}}.
#'
#' @references Mackey, L. (2009) Deflation Methods for Sparse PCA. In 
#'   \emph{Advances in Neural Information Processing Systems} (pp. 1017--1024).
#'   
#' @export
#' @param x a numeric matrix which provides the data from the first domain
#' @param xcoef a numeric data matrix with the canonical vectors related to 
#'   \code{x} as its columns.
#' @param y a numeric matrix which provides the data from the second domain
#' @param ycoef a numeric data matrix with the canonical vectors related to 
#'   \code{y} as its columns.
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
#'   
#' @return \code{acor} returns a list of class \code{nscancor} containing the 
#'   following elements: \item{cor}{the additional correlation explained by each
#'   pair of canonical variables} \item{xcoef}{copied from the
#'   input arguments}  \item{ycoef, ycenter, yscale}{copied from the input
#'   arguments} \item{xp}{the deflated data matrix corresponding to \code{x}}
#'   \item{yp}{anologous to \code{xp}}
acor <- function(x, xcoef, y, ycoef, xcenter = TRUE, ycenter = TRUE, 
                 xscale = FALSE, yscale = FALSE) {
  
  dx <- ncol(x)
  dy <- ncol(y)
  ncomp <- ncol(xcoef)
  
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
  
  W <- xcoef
  V <- ycoef
  
  corr <- rep(0, ncomp)  # additional explained correlation
  Xp <- X  # deflated X
  attr(Xp, "scaled:center") <- NULL
  attr(Xp, "scaled:scale") <- NULL
  Yp <- Y  # deflated Y
  attr(Yp, "scaled:center") <- NULL
  attr(Yp, "scaled:scale") <- NULL
  for (cc in seq_len(ncomp)) {
    w <- W[ ,cc]
    v <- V[ ,cc]
    corr[cc] <- cor(Xp%*%w, Yp%*%v)
    
    # deflate data matrices
    qx <- t(Xp)%*%(X%*%w)
    qx <- qx/normv(qx)
    Xp <- Xp - Xp%*%qx%*%t(qx)  
    
    qy <- t(Yp)%*%(Y%*%v)
    qy <- qy/normv(qy)
    Yp <- Yp - Yp%*%qy%*%t(qy)
  }
  
  nscc <- list(cor = corr, xcoef = xcoef, ycoef = ycoef,
               xcenter = if(is.null(xcen)) rep.int(0, dx) else xcen,  # return value follows cancor interface
               xscale = if(is.null(xsc)) FALSE else xsc,
               ycenter = if(is.null(ycen)) rep.int(0, dy) else ycen,
               yscale = if(is.null(ysc)) FALSE else ysc,
               xp = Xp, yp = Yp)
  class(nscc) <- "nscancor"
  return(nscc)
}