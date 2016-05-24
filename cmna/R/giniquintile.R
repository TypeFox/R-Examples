## Copyright (c) 2016, James P. Howard, II <jh@jameshoward.us>
##
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are
## met:
##
##     Redistributions of source code must retain the above copyright
##     notice, this list of conditions and the following disclaimer.
##
##     Redistributions in binary form must reproduce the above copyright
##     notice, this list of conditions and the following disclaimer in
##     the documentation and/or other materials provided with the
##     distribution.
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

#' @title Gini coefficients
#'
#' @description
#'
#' Calculate the Gini coefficient from quintile data
#'
#' @param L vector of percentages at 20th, 40th, 60th, and 80th
#' percentiles
#'
#' @details
#'
#' Calculate the Gini coefficient given the quintile data.
#'
#' @return the estimated Gini coefficient
#'
#' @family integration
#' @family newton-cotes
#'
#' @examples
#' L <- c(4.3, 9.8, 15.4, 22.7)
#' giniquintile(L)
#'
#' @references
#'
#' Leon Gerber, "A Quintile Rule for the Gini Coefficient",
#' \emph{Mathematics Magazine}, 80:2, April 2007.
#'
#' @export
giniquintile <- function(L) {
    x <- c(.2, .4, .6, .8)
    L <- x - cumsum(L / 100)

    return(25 / 144 * (3*L[1] + 2*L[2] + 2*L[3] + 3*L[4]))
}
