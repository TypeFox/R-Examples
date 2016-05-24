## Copyright (C) 2006 EPRI Solutions, Inc.
## by Tom Short, tshort@eprisolutions.com
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#for (filename in dir(pattern=".*R$")) {cat(filename);source(filename)}

filter <- function(filt, ...) UseMethod("filter")

# Octave/Matlab-compatible filter function
# y = filter (b, a, x)
filter.default <- function(filt, a, x, init, init.x, init.y, ...) {
  if(missing(init.x)) 
    init.x <- c(rep(0, length(filt) - 1))
  if(length(init.x) != length(filt) - 1)
    stop("length of init.x should match filter length-1 = ", length(filt)-1)
  if(missing(init) && !missing(init.y))
    init <- rev(init.y)
  if(all(is.na(x)))
    return(x)
  if (length(filt)) {
    x1 <- stats::filter(c(init.x, x), filt / a[1], sides = 1)
    if(all(is.na(x1)))
        return(x)
    x <- na.omit(x1, filt / a[1] , sides = 1)
    }
  if (length(a) >= 2)
    x <- stats::filter(x, -a[-1] / a[1], method = "recursive", init = init)
  x
}

filter.Arma <- function(filt, x, ...) # IIR
  filter(filt$b, filt$a, x, ...)

filter.Ma <- function(filt, x, ...) # FIR
  filter(unclass(filt), 1, x, ...)

filter.Zpg <- function(filt, x, ...) # zero-pole-gain form
  filter(as.Arma(filt), x)

MedianFilter <- function(n = 3) {
  res <- list(n = n)
  class(res) <- "MedianFilter"
  res
}

filter.MedianFilter <- function(filt, x, ...) {
  runmed(x, filt$n)
}

medfilt1 <- function(x, n = 3, ...) {
  .Deprecated("runmed", package="signal", "'medfilt1' is deprecated. Use 'runmed' of the 'stats' package instead.")
  runmed(x, n, ...)
}

spencerFilter <- function() {
  Ma(c(-3, -6, -5, 3, 21, 46, 67, 74, 67, 46, 21, 3, -5, -6, -3) / 320)
}

spencer  <- function(x)  {
    retval <- fftfilt(c(-3, -6, -5, 3, 21, 46, 67, 74, 67, 46, 21, 3, -5, -6, -3) / 320, x)
    retval <- c(rep(NA, 7), retval[-(1:14)], rep(NA, 7))
    retval
}

FilterOfOrder <- function(n, Wc, type, ...) {
  res = list(n = n, Wc = Wc, type = type, ...)
  class(res) = "FilterOfOrder"
  res
}

an <- function(degrees) {
  exp(1i*degrees*pi/180)
}

roots <- function(x, method = c("polyroot", "eigen")) {
    method <- match.arg(method)
    if(method=="polyroot")
        return(polyroot(rev(x)))
    if(!is.numeric(x))
        stop("x must be numeric")
    success <- requireNamespace("pracma", quietly = TRUE)
    if(!success)
        stop("method 'eigen' is only available if package 'pracma' is installed")
    rev(roots(as.numeric(x)))
}

Arma <- function(b, a) {
  res <- list(b = b, a = a)
  class(res) <- "Arma"
  res
}

as.Arma <- function(x, ...) UseMethod("as.Arma")

as.Arma.Arma <- function(x, ...) x

Zpg <- function(zero, pole, gain) {
  res <- list(zero = zero, pole = pole, gain = gain)
  class(res) <- "Zpg"
  res
}

as.Zpg <- function(x, ...) UseMethod("as.Zpg")

as.Zpg.Zpg <- function(x, ...) x

as.Zpg.Arma <- function(x, ...) {
  Zpg(pole = roots(x$a), zero = roots(x$b), gain = x$b[1] / x$a[1])
}

as.Zpg.Ma <- function(x, ...) {
  as.Zpg(as.Arma(x))
}

as.Arma.Zpg <- function(x, ...) {
  b = Re(x$gain*poly(x$zero))
  a = Re(poly(x$pole))
  Arma(b = b, a = a)
}

as.Arma.Ma <- function(x, ...) {
  Arma(b = unclass(x), a = 1)
}

Ma <- function(b) {
  class(b) <- "Ma"
  b
}


polyval <- function(coef, z) {
    lz <- length(z)
    if(!lz) return(numeric(0))
    n <- length(coef)
    if(!n){
        z[] <- 0
        return(z)
    }
    if(!(mode(coef) == "numeric") && !(mode(coef) == "complex"))
        stop("Argument 'coef' must be a real or complex vector.")

    ## Vectorized
    d_z <- dim(z)
    dim(z) <- lz
    y <- outer(z, (n-1):0, "^") %*% coef
    dim(y) <- d_z
    return(y)
}





conv <- function(x, y) {
  n <- length(x) + length(y) - 1
  res <- ifft(fft(postpad(x, n)) * fft(postpad(y, n)))
  if (is.numeric(x) && is.numeric(y))
    res <- Re(res)
  res
}

postpad <- function(x, n) {
  nx <- length(x)
  if (n > nx)
    c(x, rep(0,n - nx))
  else
    x[1:n]
}  

ifft <- function(x)
  fft(x, inverse = TRUE) / length(x)

sinc <- function(x){
    ifelse(x==0, 1, sin(pi*x) / (pi*x))
}

logseq <- function(from, to, n = 500)
  exp(seq(log(abs(from)), log(abs(to)), length = n))
