schweder <- function(p, xlab = "Rank of p", ylab = "p",
   drawline = NULL,
#   bh.control =list(),
   bh.lwd = 1, bh.lty = "solid", bh.col = "black",
   ls.control =list(frac = NULL),
   ls.lwd = 1, ls.lty = "dotted", ls.col = "black",
   ab.control =list(a = NULL, b = NULL),
   ab.lwd = 1, ab.lty = "dashed", ab.col = "black",
    ...) {
   bx1y1toa <- function(b, x1, y1) { # find intercept from b, x1, y1
      a <- y1 - b * x1
      a
   }
   bh <- function(p) { # parameters for Benjamini and Hochberg lowet slope line
      n <- length(p)
      rankp <- n:1
      si <- (1 - p) / (n - rankp + 1) # assumes p in non-ascending order
      di <- which(diff(si) > 0)
      ldi <- max(1, di[length(di)]) # in case no slope change
      b <- si[ldi]
      res <- list(a = bx1y1toa(b, n + 1, 1), b = b, si = si, di = di, ldi = ldi)
      res
   }
   n <- length(p)
   stopifnot(n > 1)
   toplot <- rev(sort(p))
   keep <- (toplot >= 0) & (toplot <= 1)
   toplot <- toplot[keep]
   if(length(toplot) != n) warning('Some p values omitted')
   n <- length(toplot)
   if(n < 1) stop("Not even one to plot")
   x <- n:1
   plot(x, toplot, ylab = ylab, xlab = xlab,
      ylim = c(0, 1), xlim = c(0, n+1), ...)
   res <- list(p = toplot)
   if("bh" %in% drawline) {
      params <- bh(toplot)
      abline(a = params$a, b = params$b, lty = bh.lty,
          lwd = bh.lwd, col = bh.col)
      res <- c(res, bh.params = list(a = params$a, b = params$b,
         si = params$si, di = params$di, ldi = params$ldi))
   }
   if("ls" %in% drawline) {
      if(is.null(ls.control$frac)) {
         warning("Must specify fraction for ls")
      } else {
         frac <- ls.control$frac
         newy <- toplot - 1  # transform so that line is forced
         newx <- x - (n + 1) # through (n+1),1 rather than 0,0
         plim <- newy[n * frac]
         fit <- lm(newy ~ newx - 1, subset = newy > plim)
         b <- coef(fit)[1] # correct slope
         a <- bx1y1toa(b, (n + 1), 1) # passing through (n+1,1
         abline(a = a, b = b, lty = ls.lty,
            lwd = ls.lwd, col = ls.col)
         res <- c(res, ls.params = list(a = a, b = b))
      }
   }
   if ("ab" %in% drawline) {
      if( is.null(ab.control$a) | is.null(ab.control$b) ) {
         warning("No parameters for a and b line")
      } else {
         abline(a = ab.control$a, b = ab.control$b,lty = ab.lty,
            lwd = ab.lwd, col = ab.col)
         res <- c(res, ab.params = list(a = ab.control$a, b = ab.control$b))
      }
   }
   invisible(res)
}
plot.metap <- function(x, ...) {
   schweder(x$validp, ...)
   invisible(x)
}
