#-----------------------------------------------------------------------
#     Copyright (C) 2012-2014  Serge Iovleff, University Lille 1, Inria
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as
#    published by the Free Software Foundation; either version 2 of the
#    License, or (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public
#    License along with this program; if not, write to the
#    Free Software Foundation, Inc.,
#    59 Temple Place,
#    Suite 330,
#    Boston, MA 02111-1307
#    USA
#
#    Contact : S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
#
#-----------------------------------------------------------------------
#' rtkpp is a R/STK++ bridge
#'
#' This package contains the header files for the STK++ library.
#' The typical usage is to install this package and list it in
#' the \env{LinkingTo: } line in the \file{DESCRIPTION} file of
#' other packages.
#'
#' As described at the STK++ project's home page, \url{http://www.stkpp.org},
#' STK++ is a versatile, fast, reliable and elegant collection of C++ classes
#' for statistics, clustering, linear algebra, arrays (with an Eigen-like API),
#' regression, dimension reduction, etc.
#'
#' \tabular{ll}{
#'   Package: \tab rtkpp\cr
#'   Type: \tab Package\cr
#'   Version: \tab 0.8.1\cr
#'   Date: \tab 2014-07-05\cr
#'   License: \tab GPL for the rtkpp side, LGPL for the stkpp side  + file LICENSE\cr
#'   LazyLoad: \tab yes\cr
#' }
#'
#' @rdname rtkpp-package
#' @name rtkpp
#' @aliases rtkpp
#' @docType package
#' @keywords STK++, stkpp
#' @import Rcpp
#'
#' @author
#' Author: Serge Iovleff \email{contact@@stkpp.org}
#'
#' @useDynLib rtkpp
NULL
