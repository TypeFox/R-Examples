#' Plot an optimized sample configuration
#' 
#' Plot the evolution of the energy state and the optimized sample configuration
#' 
#' @inheritParams optimACDC
#' 
#' @param x Object of class \code{OptimizedSampleConfiguration} returned by one of the 
#' \code{optim}-functions.
#' 
#' @param which Which plot should be produced: evolution of the energy state (1), optimized sample 
#' configuration (2), or both (1:2)? Defaults to \code{which = 1:2}.
#' 
#' @param boundary Object of class \code{Spatial} defining the boundary of the sampling region.
#' 
#' @param ... Other options passed to \code{plot}.
#' 
#' @rdname plot-method
#' @export
#' @method plot OptimizedSampleConfiguration
#' @aliases plot plot.OptimizedSampleConfiguration
#' @examples 
#' require(sp)
#' data(meuse.grid)
#' candi <- meuse.grid[, 1:2]
#' covars <- meuse.grid[, 5]
#' schedule <- scheduleSPSANN(initial.temperature = 5, chains = 1,
#'                            x.max = 1540, y.max = 2060, x.min = 0, 
#'                            y.min = 0, cellsize = 40)
#' set.seed(2001)
#' res <- optimCORR(points = 10, candi = candi, covars = covars, 
#'                  use.coords = TRUE, schedule = schedule)
#' plot(res)
# MAIN FUNCTION - PLOT OSC #####################################################
plot.OptimizedSampleConfiguration <-
  function (x, which = 1:2, boundary, ...) {
    
    # Do not try to plot the energy states if they have not been tracked
    if (nrow(x$objective$energy) == 2) { which <- 2 }
    
    par0 <- graphics::par()
    on.exit(suppressWarnings(graphics::par(par0)))
    if (all(which == 1:2)) { graphics::par(mfrow = c(1, 2)) }
    
    # Plot the energy states
    if (all(which == 1:2)) {
      k <- x$spsann$chains[2:3]
      k <- as.numeric(k[1] * k[2] * nrow(x$points))
      a <- x$objective$energy
      
      l <- colnames(a)
      n <- ncol(a)
      # col <- c("red", rep("black", n - 1))
      col <- c("red", grDevices::gray(seq(0, 0.5, length.out = n - 1)))
      if (n > 2) { ylim <- range(sapply(a, max)) } else { ylim <- range(a) }
      # ylim <- range(sapply(a, max))
      graphics::plot(
        1, type = 'n', xlim = c(0, k), # ylim = c(0, max(sapply(a, max)) * 1.1), 
        ylim = ylim, xlab = "jitter", ylab = "energy state")
      graphics::legend("topright", legend = l, lwd = 1, lty = 1:n, col = col)
      for(i in 1:ncol(a)) {
        graphics::lines(a[, i] ~ c(0:k), type = "l", lty = i, col = col[i])
      }
    }
    
    # Plot optimized sample configuration
    if (which == 1:2 || which == 2) {
      if (!missing(boundary)) {
        bb <- sp::bbox(boundary)
        if (methods::is(boundary, "SpatialPoints")) {
          sp::plot(x = boundary, pch = 20, cex = 0.1)
        } else {
          sp::plot(x = boundary)
        }
        graphics::points(x$points[, "x"], x$points[, "y"], pch = 20, cex = 0.5)

      } else {
        graphics::plot(x$points[, c("x", "y")], pch = 20, cex = 0.5, asp = 1)
      }
    }
  }
