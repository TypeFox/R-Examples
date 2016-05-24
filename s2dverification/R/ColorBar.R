ColorBar <- function(brks, cols = NULL, vert = TRUE, subsampleg = 1,
                     cex = 1) {
  #
  #
  #  Input arguments
  # ~~~~~~~~~~~~~~~~~
  #
  if (is.null(cols) == TRUE) {
    nlev <- length(brks) - 1
    cols <- rainbow(nlev)
  } else {
    if (length(cols) != (length(brks) - 1)) {
      stop("Inconsistent colour levels / list of colours")
    }
  }
  
  # 
  #  Plotting colorbar
  # ~~~~~~~~~~~~~~~~~~~
  #
  if (vert) {
    par(mar = c(1, 1, 1, 1.5 *( 1 + cex)), mgp = c(1, 1, 0), las = 1, cex = 1.2)
    image(1, c(1:length(cols)), t(c(1:length(cols))), axes = FALSE, col = cols, 
          xlab = '', ylab = '')
    box()
    axis(4, at = seq(0.5, length(brks) - 0.5, subsampleg), tick = TRUE, 
         labels = brks[seq(1, length(brks), subsampleg)], cex.axis = cex)
  } else {
    par(mar = c(0.5 + cex, 1, 1, 1), mgp = c(1.5, max(c(0.3,0.8*(cex-0.625))), 0),
        las = 1, cex = 1.2)
    image(1:length(cols), 1, t(t(1:length(cols))), axes = FALSE, col = cols,
          xlab = '', ylab = '')
    box()
    axis(1, at = seq(0.5, length(brks) - 0.5, subsampleg), 
         labels = brks[seq(1, length(brks), subsampleg)], cex.axis = cex)
  }
}
