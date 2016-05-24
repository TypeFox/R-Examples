##
## Copyright (c) 2010, 2011, Brandon Whitcher
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
##     * Neither the name of Rigorous Analytics Ltd. nor the names of
##       its contributors may be used to endorse or promote products 
##       derived from this software without specific prior written 
##       permission.
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
## $Id: $
##

#' Convert Decimal to Base N Number in String
#' 
#' This function converts the nonnegative integer to the specified base.
#' 
#' This function converts the nonnegative integer \code{n} to the specified
#' base, where \code{n} must be a nonnegative integer smaller than \eqn{2^52},
#' \code{base} must be an integer between 2 and 36 and \code{len} suggests the
#' length of the character string.
#' 
#' @aliases dec2base dec2hex
#' @param n Non-negative integer.
#' @param base Number between 2 and 36.
#' @param len Length of the character string.
#' @return The returned argument is a string.
#' @author Brandon Whitcher \email{bwhitcher@@gmail.com}
#' @keywords misc
#' @examples
#' 
#' x <- dec2base(23, 2)
#' 
#' @export dec2base
dec2base <- function(n, base, len=0) {
  symbols <- c(as.character(0:9), LETTERS)
  max.len <- max(trunc(log(max(n, 1)) / log(base)) + 1, len)
  ## determine digits for each number
  power <- rep(1, length(n)) * base^((max.len-1):0)
  n <- n * rep(1, max.len)
  digits <- floor((n %% (base * power)) / power)
  ## convert digits to symbols
  paste(symbols[digits + 1], collapse="")
}
#' @rdname dec2base
#' @export dec2hex
dec2hex <- function(n, len=0) {
  dec2base(n, 16, len)
}

