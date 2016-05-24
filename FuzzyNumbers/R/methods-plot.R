## This file is part of the FuzzyNumbers library.
##
## Copyright 2012-2014 Marek Gagolewski
##
##
## FuzzyNumbers is free software: you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## FuzzyNumbers is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU Lesser General Public License for more details.
##
## You should have received a copy of the GNU Lesser General Public License
## along with FuzzyNumbers. If not, see <http://www.gnu.org/licenses/>.



#' @title
#' Plot a Fuzzy Number
#'
#' @description
#' The function aims to provide a similar look-and-feel to the
#' built-in \code{\link{plot.default}} and \code{\link{curve}} function.
#'
#' @details
#' Note that if \code{from > a1} then it is set to \code{a1}.
#'
#' @param x a fuzzy number
#' @param y not used
#' @param from numeric;
#' @param to numeric;
#' @param n numeric; number of points to probe
#' @param at.alpha numeric vector; give exact alpha-cuts at which linear interpolation should be done
#' @param draw.membership.function logical; you want membership function (\code{TRUE}) or alpha-cuts plot (\code{FALSE})?
#' @param draw.alphacuts logical; defaults \code{!draw.membership.function}
#' @param xlab character; x-axis label
#' @param ylab character; y-axis label
#' @param xlim numeric;
#' @param ylim numeric;
#' @param type character; defaults \code{"l"}; plot type, e.g.~\code{"l"} for lines, \code{"p"} for points, or \code{"b"} for both
#' @param col see \code{\link{plot.default}}
#' @param lty see \code{\link{plot.default}}
#' @param pch see \code{\link{plot.default}}
#' @param lwd see \code{\link{plot.default}}
#' @param shadowdensity numeric; for shadowed sets;
#' @param shadowangle numeric; for shadowed sets;
#' @param shadowcol color specification, see \code{\link{plot.default}}; for shadowed sets;
#' @param shadowborder numeric; for shadowed sets;
#' @param add logical; add another FuzzyNumber to existing plot?
#' @param ... further arguments passed to \code{\link{plot.default}}
#'
#' @return Returns nothing really interesting.
#'
#' @exportMethod plot
#' @docType methods
#' @name plot
#' @family FuzzyNumber-method
#' @family PiecewiseLinearFuzzyNumber-method
#' @family TrapezoidalFuzzyNumber-method
#' @family DiscontinuousFuzzyNumber-method
#'
#' @rdname plot-methods
#'
#'
#' @aliases plot,FuzzyNumber,missing-method
#'          plot,TrapezoidalFuzzyNumber,missing-method
#'          plot,PiecewiseLinearFuzzyNumber,missing-method
#'          plot,DiscontinuousFuzzyNumber,missing-method
#' @usage
#' \S4method{plot}{FuzzyNumber,missing}(x, y, from=NULL, to=NULL, n=101, at.alpha=NULL,
#' draw.membership.function=TRUE, draw.alphacuts=!draw.membership.function,
#' xlab=NULL, ylab=NULL, xlim=NULL, ylim=NULL,
#' type="l", col=1, lty=1, pch=1, lwd=1,
#' shadowdensity=15, shadowangle=45, shadowcol=col, shadowborder=NULL,
#' add=FALSE, ...)
#'
#' \S4method{plot}{TrapezoidalFuzzyNumber,missing}(x, y, from=NULL, to=NULL,
#' draw.membership.function=TRUE, draw.alphacuts=!draw.membership.function,
#' xlab=NULL, ylab=NULL, xlim=NULL, ylim=NULL,
#' type="l", col=1, lty=1, pch=1, lwd=1, add=FALSE, ...)
#'
#' \S4method{plot}{PiecewiseLinearFuzzyNumber,missing}(x, y, from=NULL, to=NULL,
#' draw.membership.function=TRUE, draw.alphacuts=!draw.membership.function,
#' xlab=NULL, ylab=NULL, xlim=NULL, ylim=NULL,
#' type="l", col=1, lty=1, pch=1, lwd=1, add=FALSE, ...)
#'
#' \S4method{plot}{DiscontinuousFuzzyNumber,missing}(x, y, from=NULL, to=NULL,
#' n=101, draw.membership.function=TRUE, draw.alphacuts=!draw.membership.function,
#' xlab=NULL, ylab=NULL, xlim=NULL, ylim=NULL,
#' type="l", col=1, lty=1, pch=1, lwd=1,
#' add=FALSE, ...)
#'
#' @examples
#' plot(FuzzyNumber(0,1,2,3), col="gray")
#' plot(FuzzyNumber(0,1,2,3, left=function(x) x^2, right=function(x) 1-x^3), add=TRUE)
#' plot(FuzzyNumber(0,1,2,3, lower=function(x) x, upper=function(x) 1-x), add=TRUE, col=2)
if (!isGeneric("plot"))
   setGeneric("plot", function(x, y, ...) standardGeneric("plot"))



setMethod(
   f="plot",
   signature(x="FuzzyNumber", y="missing"),
   definition=function(x, y, from=NULL, to=NULL, n=101, at.alpha=NULL,
      draw.membership.function=TRUE, draw.alphacuts=!draw.membership.function,
      xlab=NULL, ylab=NULL, xlim=NULL, ylim=NULL,
      type="l", col=1, lty=1, pch=1, lwd=1,
      shadowdensity=15, shadowangle=45, shadowcol=col, shadowborder=NULL,
      add=FALSE, ...)
   {
      draw.membership.function <- identical(draw.membership.function, TRUE);
      draw.alphacuts <- identical(draw.alphacuts, TRUE);

      if (!draw.membership.function && !draw.alphacuts)
         stop("Provide `draw.alphacuts' or `draw.membership.function'");

      if (!draw.alphacuts && (is.null(ylim) || !is.numeric(ylim) || length(ylim) != 2))
         ylim <- c(0,1);

      if ( draw.alphacuts && (is.null(xlim) || !is.numeric(xlim) || length(xlim) != 2))
         xlim <- c(0,1);

      if (!draw.alphacuts && is.null(xlab))
         xlab <- "x";

      if (!draw.alphacuts && is.null(ylab))
         ylab <- expression(alpha);

      if ( draw.alphacuts && is.null(xlab))
         xlab <- expression(alpha);

      if ( draw.alphacuts && is.null(ylab))
         ylab <- "x";

      drawX     <- !(is.na(x@left(0)));
      drawAlpha <- !(is.na(x@lower(0)));

      add <- identical(add, TRUE);

      if (dev.cur() == 1L && add)
      {
         warning("`add' will be ignored as there is no existing plot")
         add <- FALSE;
      }

      if (n <= 0) n <- 0;

      if (draw.alphacuts)
      {
## -----------------------------------------------------------------------
## ----------------------------------------------  draw.alphacuts --------

         # prepare ylim, from, to
         if (is.null(from) || is.null(to))
         {
            ylim <-
               if (is.numeric(ylim) && length(ylim) == 2)
               {
                  ylim;
               } else if (add)
               {
                  usr <- par("usr")[3L:4L]
                  if (par("yaxs") == "r") usr <- extendrange(usr, f = -1.0/27.0);
                  usr
               } else
               {
                  extendrange(c(x@a1, x@a4), f = 1.0/27.0); # core + epsilon
               }

            if (is.null(from)) from <- ylim[1L];
            if (is.null(to))   to <- ylim[2L];
         } else if (is.null(ylim))
         {
            ylim <- c(from, to);
         }

         if (from > x@a1) from <- x@a1;
         if (to   < x@a4) to   <- x@a4;


         if (!drawX && !drawAlpha)
         {
## ========================================= draw.alpha: shadowed FN
            matplot(NA, NA, type=type,
               xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim,
               col=col, lty=lty, pch=pch, lwd=lwd, add=add, ...);
            rect(0, x@a1, 1, x@a2, density=shadowdensity,
               col=shadowcol, angle=shadowangle, border=shadowborder);
            rect(1, x@a3, 0, x@a4, density=shadowdensity,
               col=shadowcol, angle=shadowangle, border=shadowborder);
## ==================================================================
         } else
         {
            if (drawAlpha)
            {
## ========================================= draw.alpha: from alpha cuts
               if (!is.numeric(at.alpha) || is.unsorted(at.alpha) ||
                     any(at.alpha <= 0 | at.alpha >= 1) || length(at.alpha) == 0)
               {
                  at.alpha <- seq(1/(n+1),1-1/(n+1),length.out=n);
               }

               if (length(at.alpha) == 0)
               {
                  xvals1 <- numeric(0);
                  xvals2 <- numeric(0);
               } else if (length(at.alpha) == 1)
               {
                  xvals <- alphacut(x, at.alpha);
                  xvals1 <- xvals[1];
                  xvals2 <- rev(xvals[2]);
               } else
               {
                  xvals <- alphacut(x, at.alpha);
                  xvals1 <- xvals[,1];
                  xvals2 <- rev(xvals[,2]);
               }
               alpha1 <- at.alpha;
               alpha2 <- rev(at.alpha);
## =====================================================================
            } else
            {
## ================================= draw.alpha: from sides
               if (n == 0)
               {
                  xvals1 <- numeric(0);
                  xvals2 <- numeric(0);
               } else
               {
                  xvals1 <- seq(x@a1, x@a2, length.out=n+2); xvals1 <- xvals1[-c(1,n+2)];
                  xvals2 <- seq(x@a3, x@a4, length.out=n+2); xvals2 <- xvals2[-c(1,n+2)];
               }
               alpha1 <- evaluate(x, xvals1);
               alpha2 <- evaluate(x, xvals2);
## =========================================================
            }

            draw.y1 <- c(x@a1, xvals1, x@a2);
            draw.x1 <- c(0,    alpha1, 1);
            draw.y2 <- c(x@a3, xvals2, x@a4);
            draw.x2 <- c(1,    alpha2, 0);

            # draw from sides or alpha-cuts
            matplot(cbind(draw.x1, draw.x2), cbind(draw.y1, draw.y2),
                    type=type, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, col=col,
                    lty=lty, pch=pch, lwd=lwd, add=add, ...);

## ---------------------------------------------  /draw.alphacuts --------
## -----------------------------------------------------------------------
         }
      } else if (draw.membership.function)
      {
## -----------------------------------------------------------------------
## ------------------------------------  draw.membership.function --------

         # prepare xlim, from, to
         if (is.null(from) || is.null(to))
         {
            xlim <-
               if (is.numeric(xlim) && length(xlim) == 2)
               {
                  xlim;
               } else if (add)
               {
                  usr <- par("usr")[1L:2L]
                  if (par("xaxs") == "r") usr <- extendrange(usr, f = -1.0/27.0);
                  usr
               } else
               {
                  extendrange(c(x@a1, x@a4), f = 1.0/27.0); # core + epsilon
               }

            if (is.null(from)) from <- xlim[1L];
            if (is.null(to))   to <- xlim[2L];
         } else if (is.null(xlim))
         {
            xlim <- c(from, to);
         }

         if (from > x@a1) from <- x@a1;
         if (to   < x@a4) to   <- x@a4;



         if (!drawX && !drawAlpha)
         {
## ================================= draw.membership.funtion: shadowed FN
            matplot(c(from, x@a1), c(0,0), type=type,
               xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim,
               col=col, lty=lty, pch=pch, lwd=lwd, add=add, ...);
            rect(x@a1, 0, x@a2, 1, density=shadowdensity,
               col=shadowcol, angle=shadowangle, border=shadowborder);
            rect(x@a3, 1, x@a4, 0, density=shadowdensity,
               col=shadowcol, angle=shadowangle, border=shadowborder);
            matplot(c(x@a4, to),   c(0,0), type=type,
               col=col, lty=lty, pch=pch, lwd=lwd, add=TRUE, ...);
            matplot(c(x@a2, x@a3), c(1,1), type=type,
               col=col, lty=lty, pch=pch, lwd=lwd, add=TRUE, ...);
## =======================================================================
         } else
         {
            if (drawAlpha && (!drawX || !is.null(at.alpha)))
            {
## ================================= draw.membership.funtion: from alpha cuts
               if (!is.numeric(at.alpha) || is.unsorted(at.alpha) ||
                     any(at.alpha <= 0 | at.alpha >= 1) || length(at.alpha) == 0)
               {
                  at.alpha <- seq(1/(n+1),1-1/(n+1),length.out=n);
               }

               if (length(at.alpha) == 0)
               {
                  xvals1 <- numeric(0);
                  xvals2 <- numeric(0);
                  alpha1 <- numeric(0);
                  alpha2 <- numeric(0);
               } else if (length(at.alpha) == 1)
               {
                  xvals <- alphacut(x, at.alpha);
                  xvals1 <- xvals[1];
                  xvals2 <- rev(xvals[2]);
                  alpha1 <- at.alpha;
                  alpha2 <- at.alpha;
               } else
               {
                  xvals <- alphacut(x, at.alpha);

                  xvals1 <- xvals[,1];
                  xvals2 <- rev(xvals[,2]);
                  alpha1 <- at.alpha;
                  alpha2 <- rev(at.alpha);
               }
## ===========================================================================
            } else
            {
## ================================= draw.membership.funtion: from sides
               if (n == 0)
               {
                  xvals1 <- numeric(0);
                  xvals2 <- numeric(0);
                  alpha1 <- numeric(0);
                  alpha2 <- numeric(0);
               } else
               {
                  xvals1 <- seq(x@a1, x@a2, length.out=n+2); xvals1 <- xvals1[-c(1,n+2)];
                  xvals2 <- seq(x@a3, x@a4, length.out=n+2); xvals2 <- xvals2[-c(1,n+2)];
                  alpha1 <- evaluate(x, xvals1);
                  alpha2 <- evaluate(x, xvals2);
               }
## =====================================================================
            }

            draw.x <- c(from, x@a1, xvals1, x@a2, x@a3, xvals2, x@a4, to);
            draw.y <- c(0,    0,    alpha1, 1,    1,    alpha2, 0,    0);

            # draw from sides or alpha-cuts
            matplot(draw.x, draw.y,
                    type=type, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, col=col,
                    lty=lty, pch=pch, lwd=lwd, add=add, ...);

## ------------------------------------ /draw.membership.function --------
## -----------------------------------------------------------------------
         }
      }

      invisible(NULL)
   }
)



setMethod(
   f="plot",
   signature(x="DiscontinuousFuzzyNumber", y="missing"),
   definition=function(x, y, from=NULL, to=NULL, n=101,
                       draw.membership.function=TRUE, draw.alphacuts=!draw.membership.function,
                       xlab=NULL, ylab=NULL, xlim=NULL, ylim=NULL,
                       type="l", col=1, lty=1, pch=1, lwd=1,
                       add=FALSE, ...)
   {
      at.alpha <- unique(sort(pmax(0,pmin(1,c(x@discontinuities.lower-1e-16,
                             x@discontinuities.lower+1e-16,
                             x@discontinuities.upper-1e-16,
                             x@discontinuities.upper+1e-16,
                             seq(0, 1, length.out=n))))))
      callNextMethod(x, at.alpha=at.alpha,
                     from=from, to=to, type=type, xlab=xlab, ylab=ylab,
                     xlim=xlim, ylim=ylim, col=col, lty=lty, pch=pch, lwd=lwd,
                     draw.membership.function=draw.membership.function,
                     draw.alphacuts=draw.alphacuts, add=add, ...)
   }
)




setMethod(
   f="plot",
   signature(x="TrapezoidalFuzzyNumber", y="missing"),
   definition=function(x, y, from=NULL, to=NULL,
                       draw.membership.function=TRUE, draw.alphacuts=!draw.membership.function,
                       xlab=NULL, ylab=NULL, xlim=NULL, ylim=NULL,
                       type="l", col=1, lty=1, pch=1, lwd=1,
                       add=FALSE, ...)
   {
      callNextMethod(x, n=0, at.alpha=numeric(0),
                     from=from, to=to, type=type, xlab=xlab, ylab=ylab,
                     xlim=xlim, ylim=ylim, col=col, lty=lty, pch=pch, lwd=lwd,
                     draw.membership.function=draw.membership.function,
                     draw.alphacuts=draw.alphacuts, add=add, ...)
   }
)




setMethod(
   f="plot",
   signature(x="PiecewiseLinearFuzzyNumber", y="missing"),
   definition=function(x, y, from=NULL, to=NULL,
                       draw.membership.function=TRUE, draw.alphacuts=!draw.membership.function,
                       xlab=NULL, ylab=NULL, xlim=NULL, ylim=NULL,
                       type="l", col=1, lty=1, pch=1, lwd=1,
                       add=FALSE, ...)
   {
      callNextMethod(x, at.alpha=x@knot.alpha,
                     from=from, to=to, type=type, xlab=xlab, ylab=ylab,
                     xlim=xlim, ylim=ylim, col=col, lty=lty, pch=pch, lwd=lwd,
                     draw.membership.function=draw.membership.function,
                     draw.alphacuts=draw.alphacuts, add=add, ...)
   }
)
