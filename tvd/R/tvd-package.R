# tvd-package.R: Roxygen2 documentation file for the tvd package.
# 
# Copyright 2014 Mark Pinese
#
# This file is distributed under the terms of the Eclipse Public 
# License v1.0, available at:
# https://www.eclipse.org/org/documents/epl-v10.html


#' Total Variation Denoising
#'
#' tvd is a package for Total Variation Denoising, a regularized 
#' procedure for removing noise from piecewise constant signals
#' whilst retaining edges.  Currently implements Condat's algorithm
#' for fast 1D TVD, in function tvd1d.
#'
#' @useDynLib tvd
#' @name tvd-package
#' @docType package
#' @title Total Variation Denoising
#' @author Mark Pinese \email{m.pinese@@garvan.org.au}
#' @references Condat, L. (2013) A Direct Algorithm for 1-D Total Variation Denoising.
#'   IEEE Signal Processing Letters 20(11): 1054-1057.  \url{doi:10.1109/LSP.2013.2278339}
#' @keywords package
#' @seealso \code{\link{tvd1d}}
NULL
