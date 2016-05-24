## This file is part of the 'agop' library.
##
## Copyright 2013 Marek Gagolewski, Anna Cena
##
## Parts of the code are taken from the 'CITAN' R package by Marek Gagolewski
##
## 'agop' is free software: you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## 'agop' is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU Lesser General Public License for more details.
##
## You should have received a copy of the GNU Lesser General Public License
## along with 'agop'. If not, see <http://www.gnu.org/licenses/>.


#' @title
#' Draws a Graphical Representation of a Given Vector
#'
#' @description
#' Draws a step function that represents given numeric vector
#' with elements in \eqn{[0,\infty]}.
#' 
#' @details
#' In \pkg{agop}, a given vector \eqn{x=(x_1,\dots,x_n)} can be represented by a
#' step function defined for \eqn{0\le y<n} and given by:
#' \deqn{\pi(y)=x_{(n-\lfloor y+1\rfloor+1)}}{\pi(y)=x_{(n-floor(y+1)+1)}}
#' (for \code{type == 'right.continuous'})
#' or for \eqn{0< y\le n} \deqn{\pi(y)=x_{(n-\lfloor y\rfloor+1)}}{\pi(y)=x_{(n-floor(y)+1)}}
#' (for \code{type == 'left.continuous'}, the default)
#' or by a curve joining the points \eqn{(0, x_{(n)})},
#' \eqn{(1, x_{(n)})}, \eqn{(1, x_{(n-1)})}, \eqn{(2, x_{(n-1)})},
#' ..., \eqn{(n, x_{(1)})}.
#' Here, \eqn{x_{(i)}} denotes the
#' \eqn{i}-th smallest value in \eqn{x}.
#' 
#' In bibliometrics, a step function of one of the two above-presented types
#' is called a citation function.
#' 
#' For historical reasons, this function is also available via its alias,
#' \code{plot.citfun} [but its usage is deprecated].
#'
#' @param x non-negative numeric vector
#' @param type character; type of the graphical \code{'left.continuous'} (the default)
#' or \code{'right.continuous'} for step functions and \code{'curve'} for
#' a continuous step curve
#' @param extend logical; should the plot be extended infinitely to the right?
#' Defaults to \code{FALSE}
#' @param add logical; indicates whether to start a new plot, \code{FALSE} by default
#' @param pch,col,lty,lwd,cex,xmarg graphical parameters
#' @param col.steps,lty.steps,lwd.steps graphical parameters, used only
#' for \code{type} of \code{'left.continuous'} and \code{'right.continuous'} only
#' @param ylim,xlim,xlab,ylab,main,... additional graphical parameters,
#' see \code{\link{plot.default}}
#' @return
#' nothing interesting
#' 
#' @examples
#' john_s <- c(11,5,4,4,3,2,2,2,2,2,1,1,1,0,0,0,0)
#' plot_producer(john_s, main="Smith, John", col="red")
#' @export
#' @name plot_producer
#' @rdname plot_producer
#' @family visualization
#' @aliases plot.citfun
plot_producer <- function(x,
   type=c('left.continuous', 'right.continuous', 'curve'),
   extend=FALSE, add=FALSE,
   pch=1, col=1, lty=1, lwd=1, cex=1,
   col.steps=col, lty.steps=2, lwd.steps=1,
   xlab="", ylab="", main="",
   xmarg=10, xlim=c(0,length(x)*1.2), ylim=c(0, max(x)), 
   ...)
{
   stopifnot(length(x) > 0, is.numeric(x), all(x >= 0))
   n <- length(x)
   x <- sort(x, decreasing=T)
   
   type <- match.arg(type)
   extend <- identical(extend, TRUE)
   add <- identical(add, TRUE)
   
   if (!add) {
      # start a new plot
      plot(NA, NA, ylim=ylim, xlim=xlim,
         xlab=xlab, ylab=ylab, main=main, ...)
   }
   
   
   # look for unique values
   wdx <- which(diff(x) != 0) 
   px <- c(wdx, n)
   py <- c(x[wdx], x[n])
   
   if (extend) {
      if (py[length(py)] == 0) {
         px[length(px)] <- par("usr")[2]+1 # +1 == magic constant :-)
      }
      else {
         px <- c(px, par("usr")[2]+1)
         py <- c(py, 0)
      }
   }
   
   segments(c(0, px[-length(px)]), py, px, py,
      col=col, pch=pch, lty=lty, lwd=lwd, cex=cex, ...)

   if (type == 'right.continuous') {
      points(px, py, 
         col=col, pch=pch, lty=lty, lwd=lwd, cex=cex, ...)
   }
   else if (type == 'left.continuous') {
      points(c(0, px[-length(px)]), py, 
         col=col, pch=pch, lty=lty, lwd=lwd, cex=cex, ...)
   }
   
   if (type == 'curve') {
      segments(px[-length(px)], py[-length(px)], px[-length(px)], py[-1],
         col=col, lty=lty, lwd=lwd, ...)
   }
   else {
      segments(px[-length(px)], py[-length(px)], px[-length(px)], py[-1],
         col=col.steps, lty=lty.steps, lwd=lwd.steps, ...)
   }
}


#' @export
plot.citfun <- plot_producer # deprecated













# #' The \eqn{r_p}-curve appears in the definition of the \eqn{r_p}-index
# #' (see Gagolewski, Grzegorzewski, 2009) and the \code{\link{index.rp}} function.
# #'
# #' @title Draw the r_p-curve of given radius
# #' @param r radius of the \eqn{r_p}-curve; \eqn{r>0}.
# #' @param p index order, \eqn{p \in [1,\infty]}{p in [1,\infty]}; defaults \eqn{\infty} (\code{Inf}).
# #' @param n integer; the maximal number of values at which to evaluate the underlying function.
# #' @param ... additional graphical parameters.
# #' @seealso \code{\link{index.rp}}, \code{\link{plot.citfun}}, \code{\link{plot.default}}
# #' @examples
# #' john_s <- c(11,5,4,4,3,2,2,2,2,2,1,1,1,0,0,0,0);
# #' plot.citfun(john_s, main="Smith, John");
# #' curve.add.rp(index.rp(john_s), col="green");
# #' curve.add.rp(index.rp(john_s,1), 1, col="blue");
# #' curve.add.rp(index.rp(john_s,2), 2, col="red");
# #' @references
# #' Gagolewski M., Grzegorzewski P., A geometric approach to the construction of scientific impact indices, Scientometrics, 81(3), 2009a, 617-634.\cr
# #' @export
# curve.add.rp <- function(r, p=Inf, n=101, ...)
# {
#    if (length(p) != 1 || mode(p) != "numeric") stop("p must be a single numeric value");
#    if (p < 1) stop("p must be >= 1");
# 
#    if (length(r) != 1 || mode(r) != "numeric") stop("r must be a single numeric value");
#    if (r <= 0) stop("r must be > 0");
# 
#    if (is.finite(p))
#    {
#       px <- seq(0,r,length=n);
#       py <- (r^p-px^p)^(1.0/p);
#       lines(px, py, ...);
#    } else
#    {
#       lines(c(0,r), c(r,r), ...);
#       points(r, r, ...);
#    }
# }




# #' The \eqn{l_p}-curve appears in the definition of the \eqn{l_p}-index
# #' (see Gagolewski, Grzegorzewski, 2009) and the \code{\link{index.lp}} function.
# #'
# #' @title Draw the l_p-curve of given size
# #' @param ab size of the \eqn{l_p}-curve; positive numeric vector of length 2.
# #' @param p index order, \eqn{p \in [1,\infty]}{p in [1,\infty]}; defaults \eqn{\infty} (\code{Inf}).
# #' @param n integer; the maximal number of values at which to evaluate the underlying function.
# #' @param ... additional graphical parameters.
# #' @seealso \code{\link{index.lp}}, \code{\link{plot.citfun}}, \code{\link{plot.default}}
# #' @examples
# #' john_s <- c(11,5,4,4,3,2,2,2,2,2,1,1,1,0,0,0,0);
# #' plot.citfun(john_s, main="Smith, John");
# #' curve.add.lp(index.lp(john_s), col="green");
# #' curve.add.lp(index.lp(john_s,1), 1, col="blue");
# #' curve.add.lp(index.lp(john_s,2), 2, col="red");
# #' @references
# #' Gagolewski M., Grzegorzewski P., A geometric approach to the construction of scientific impact indices, Scientometrics, 81(3), 2009a, 617-634.\cr
# #' @export
# curve.add.lp <- function(ab, p=Inf, n=101, ...)
# {
#    if (length(p) != 1 || mode(p) != "numeric") stop("p must be a single numeric value");
#    if (p < 1) stop("p must be >= 1");
# 
#    if (length(ab) != 2 || mode(ab) != "numeric") stop("ab must be a numeric vector of length 2");
#    if (any(ab <= 0)) stop("ab must be > 0");
# 
#    a <- ab[1];
#    b <- ab[2];
# 
#    if (is.finite(p))
#    {
#       px <- seq(0,a,length=n);
#       py <- (b^p-(b/a*px)^p)^(1.0/p);
#       lines(px, py, ...);
#    } else
#    {
#       lines(c(0,a), c(b,b), ...);
#       points(a, b, ...);
#    }
# }
