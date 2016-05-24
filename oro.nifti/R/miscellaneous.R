##
##
## Copyright (c) 2009, Brandon Whitcher and Volker Schmid
## All rights reserved.
## 
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are
## met:
## 
##     * Redistributions of source code must retain the above copyright
##       notice, this list of conditions and the following disclaimer. 
##     * Redistributions in binary form must reproduce the above
##       copyright notice, this list of conditions and the following
##       disclaimer in the documentation and/or other materials provided
##       with the distribution.
##     * The names of the authors may not be used to endorse or promote
##       products derived from this software without specific prior
##       written permission.
## 
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
## "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
## LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
## A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
## HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
## SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
## LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
## DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
## THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
## (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
## OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
## 
##

#' Dimension Accessor Functions
#' 
#' Functions to extract the higher dimensions from ANALYZE/NIfTI data.
#' 
#' Simple calls to \code{dim} to replicate the functionality of \code{nrow} and
#' \code{ncol} for higher dimensions of an array that are commonly required
#' when manipulating medical imaging data.
#' 
#' @aliases nsli NSLI ntim NTIM
#' @param x is a three- or four-dimensional array (e.g., read in from an
#' ANALYZE/NIfTI file).
#' @return Third (slice) or fourth (time) dimension of the array.
#' @author Brandon Whitcher \email{bwhitcher@@gmail.com}
#' @seealso \code{\link{readNIfTI}}, \code{\link{readANALYZE}}
#' @rdname miscellaneous
#' @keywords misc
#' @export
nsli <- function(x) {
  dim(x)[3]
}
#' @rdname miscellaneous
#' @export
NSLI <- function(x) {
  dim(x)[3]
}
#' @rdname miscellaneous
#' @export
ntim <- function(x) {
  dim(x)[4]
}
#' @rdname miscellaneous
#' @export
NTIM <- function(x) {
  dim(x)[4]
}
