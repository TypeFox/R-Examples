## This file is part of the R package softclassval.
##
## softclassval is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## Foobar is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
##
## In R, you can display the license with
## > RShowDoc ("GPL-3")
##
## In addition to the general terms of the GPL you are required to cite
## the package appropriately (according to subsec).
## R's citation command gives the correct citation:
## > citation ("softclassval")
##
## Copyright 2010 - 2011 C. Beleites


##' Mark operator as deviation measure
##'
##' The operators measure either a performance (i.e. accordance between reference and prediction) or
##' a deviation. \code{dev (op) == TRUE} marks operators measuring deviation.
##' 
##' @param op the operator (function)
##' @return logical indicating the type of operator. \code{NULL} if the attribute is missing.
##' @author Claudia Beleites
##' @seealso \code{\link{sens}} \code{\link{post}}
##' @export 
##' @include softclassval.R
##' @include unittestdata.R
##'
##' @examples
##'
##' dev (wRMSE)
##' myop <- function (r, p) p * (r == 1)
##' dev (myop) <- TRUE
##' 

dev <- function (op)
  attr (op, "dev")

##' @usage dev (op) <- value
##' @rdname dev
##' @param value logical indicating the operator type
##' @export "dev<-"
"dev<-" <- function (op, value){
  stopifnot (is.logical (value), !is.na (value))

  attr (op, "dev") <- value

  op
}

.test (dev) <- function (){
  myop <- function (){}
  checkTrue (is.null (dev (myop)))
  dev (myop) <- TRUE
  checkTrue (dev (myop))
  dev (myop) <- FALSE
  checkTrue (!dev (myop))
  checkException (dev (myop) <- NULL)
  checkException (dev (myop) <- NA)
}
