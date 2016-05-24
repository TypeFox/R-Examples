#' @title Utilities for LC class
#' @name LC-utils
#' @aliases plot.LC summary.LC LC_coordinates2control_settings
#' @description 
#' 
#' The \code{"LC"} class is the core property of the LICORS model as it specifies
#' the spatio-temporal neighborhood of the past and future light cone. 
#' The function \code{\link{setup_LC_geometry}} generates an instance of 
#' the \code{"LC"} class.
#' 
#' \code{plot.LC} plots LCs of \eqn{(1+1)D} and \eqn{(2+1)D} systems with 
#' respect to the origin \eqn{(\mathbf{r}, t) = (\boldsymbol 0, 0)}. 
#' This is especially useful for a quick check if the LC geometry 
#' specified by \code{\link{setup_LC_geometry}} is really the intended one.
#'
#' \code{summary.LC} prints a summary of the LC geometry.  
#' Returns (invisible) the summary matrix.
#'
#' \code{LC_coordinates2control_setting} computes auxiliary measures given the LC 
#' geometry such as horizon, spatial/temporal extension, etc. This function 
#' should not be called by the user directly; only by \code{\link{get_spacetime_grid}}.
#' 
NULL

#' @rdname LC-utils
#' @keywords hplot
#' @method plot LC
#' @param x an object of class \code{"LC"} (see \code{\link{setup_LC_geometry}})
#' @param cex.axis The magnification to be used for axis annotation relative to 
#' the current setting of \code{cex}.
#' @param cex.lab The magnification to be used for x and y labels 
#' relative to the current setting of \code{cex}.
#' @param ... optional arguments passed to \code{plot}.
#' @export
#' @examples
#' 
#' aa = setup_LC_geometry(horizon = list(PLC = 2, FLC= 1), speed=1, space.dim = 1, shape = "cone")
#' plot(aa)
#' bb = setup_LC_geometry(horizon = list(PLC = 2, FLC= 1), speed=1, space.dim = 1, shape = "revcone")
#' plot(bb)
#' 

plot.LC <- function(x, cex.axis = 2, cex.lab = 2, ...) {
  object <- x
  
  LC.coordinates <- rbind(object$coordinates$PLC, object$coordinates$FLC)
  if (object$space.dim == 0) {
    plot(cbind(LC.coordinates, 0), type = "n", cex.lab = cex.lab, xlab = "Time", 
         axes = FALSE, xlim = range(LC.coordinates) + c(-0.5, 0.5), ylab = "",
         ylim = c(-1, 1), ...)
    box()
    abline(v = 0)
    axis(1, at = unique(LC.coordinates[, 1]), cex.axis = cex.axis)
    grid()
    points(cbind(object$coordinates$FLC, 0), pch = "+", cex = 3, col = "blue")
    points(cbind(object$coordinates$PLC, 0), pch = "--", cex = 3, col = "red")
    text(mean(object$coordinates$PLC), -0.5, "PLC", col = "red", cex = 2)
    text(mean(object$coordinates$FLC), 0.5, "FLC", col = "blue", cex = 2)
  } else if (object$space.dim == 1) {
    plot(LC.coordinates, type = "n", cex.lab = cex.lab, xlab = "Time", ylab = "Space", 
         axes = FALSE,  xlim = range(LC.coordinates[, 1]) + c(-0.5, 0.5), ...)
    box()
    abline(v = 0)
    
    axis(1, at = unique(LC.coordinates[, 1]), cex.axis = cex.axis)
    axis(2, at = unique(LC.coordinates[, 2]), cex.axis = cex.axis)
    grid()
    points(object$coordinates$FLC, pch = "+", cex = 3, col = "blue")
    points(object$coordinates$PLC, pch = "--", cex = 3, col = "red")
    text(max(object$coordinates$PLC[, 1]) + 0.25, max(LC.coordinates[, 2]) - 0.5, "PLC", col = "red", cex = 2)
    text(min(object$coordinates$FLC[, 1]) + 0.25, min(LC.coordinates[, 2]) + 0.5, "FLC", col = "blue", cex = 2)
  } else if (object$space.dim == 2) {
    print("Not yet implemented.")
  } else {
    plot.new()
    legend("top", paste("Can't plot (3+1)D systems.", cex = 2))
  }
} 

#' @rdname LC-utils
#' @keywords models print
#' @param object an object of class \code{"LC"}
#' @param verbose logical; if \code{TRUE} LC information is printed in the 
#' console
#' @export
#' @method summary LC
#' @examples
#' 
#' aa = setup_LC_geometry(horizon = list(PLC = 2, FLC = 0), speed=1, space.dim = 1, shape = "cone")
#' summary(aa)
#' 
summary.LC <- function(object, verbose = TRUE, ...) {
  
  tmp.mat <- rbind(object$horizon)
  tmp.mat <- rbind(tmp.mat, object$speed)
  tmp.mat <- as.data.frame(tmp.mat)
  tmp.mat <- rbind(tmp.mat, object$shape)
  tmp.mat <- rbind(tmp.mat, c(object$n.p, object$n.f))
  
  print(tmp.mat)
  colnames(tmp.mat) <- c("PLC", "FLC")
  rownames(tmp.mat) <- c("horizon", "speed", "shape", "dimensionality")
  
  if (verbose) {
    cat(rep("*", 20))
    cat("\n \n")
    cat(paste("The field extends over a ", object$space.dim, "-dimensional space. \n", 
              sep = ""))
    cat("Light cones have therefore the following characteristics: \n \n")
    print(tmp.mat)
    cat("\n")
    cat(rep("*", 20))
  }
  
  out <- list(table = tmp.mat)
  class(out) <- "summary.LC"
  invisible(out)
}

#' @rdname LC-utils
#' @param LC.coordinates template of a light cone (with respect to origin)
#' @keywords method
#' @export
#' @seealso \code{\link{compute_LC_coordinates}}
#' @examples
#' 
#' aa = setup_LC_geometry(horizon = list(PLC = 2, FLC = 0), speed=1, space.dim = 1, shape = "cone")
#' LC_coordinates2control_settings(aa$coordinates)
#' 

LC_coordinates2control_settings <- function(LC.coordinates) {
  
  out <- list(n.p = nrow(LC.coordinates$PLC),
              n.f = nrow(LC.coordinates$FLC), 
              space.dim = ncol(LC.coordinates$FLC) - 1)
  out$horizon <- list(PLC = max(-LC.coordinates$PLC[, 1]),
                      FLC = max(LC.coordinates$FLC[, 1]))

  joint.LC.coordinates <- rbind(LC.coordinates$PLC, LC.coordinates$FLC)
  if (out$space.dim == 0) {
    out$space.cutoff <- 0
  } else if (out$space.dim == 1) {
    out$space.cutoff <- max(joint.LC.coordinates[, -1])
  } else {
    out$space.cutoff <- apply(joint.LC.coordinates[, -1], 2, max)
  }
  return(out)
}

