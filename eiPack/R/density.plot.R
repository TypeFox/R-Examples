density.plot <- function(x, by = "column", col,
                         xlim = c(0,1), ylim,
                         main = "", sub = NULL, xlab,
                         ylab, lty = par("lty"), 
                         lwd = par("lwd"), ...) {
  UseMethod("density.plot", x)
}
