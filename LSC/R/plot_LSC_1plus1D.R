#' @rdname LSC-utils
#' @description
#' \code{\link{plot_LSC_1plus1D}} plots LSC for a (1+1)D field.
#' @keywords hplot
#' @param widths passed to \code{\link[graphics]{layout}} for dividing the plotting
#' region horizontally. A vector of length 2: image (left) & temporal 
#' average (right)
#' @export
#' @examples
#' \dontrun{
#' data(contCA00)
#' 
#' temp_lsc = states2LSC(states = contCA00$predictive_states - min(contCA00$predictive_states) + 1)
#' class(temp_lsc) = c("LSC", "LSC_1plus1D")
#' plot_LSC_1plus1D(temp_lsc)
#' }
#'

plot_LSC_1plus1D = function(z, col = NULL, lsc.unit = "bits", 
                            heights = c(2,5), widths = c(5,2)) {
  
  LSC <- z

  if (is.null(col)) {
    col <- colorRampPalette(brewer.pal(9, name="YlOrRd"))(100)
  } 
  if (length(col) == 1 & is.character(col)){
    col <- colorRampPalette(brewer.pal(9, name=col))(100)
  }
  
  LSC_temporal_avg <- colMeans(LSC)
  LSC_spatial_avg <- rowMeans(LSC)
  TT <- nrow(LSC)
  card_space <- ncol(LSC)
  
  #LSC_temporal_avg_smooth <- gam(LSC_temporal_avg ~ s(c(1:card_space)))$fitted.values
  #LSC_spatial_avg_smooth <- gam(LSC_spatial_avg ~ s(c(1:TT)))$fitted.values
  
  LSC_temporal_avg_smooth <- loess(LSC_temporal_avg ~ c(1:card_space), span = 0.2)$fit
  LSC_spatial_avg_smooth <- loess(LSC_spatial_avg ~ c(1:TT), span = 0.2)$fit 
  
  #LSC_temporal_avg_smooth <- loess.as(c(1:card_space), LSC_temporal_avg, degree = 2)$fit
  #LSC_spatial_avg_smooth <- loess.as(c(1:TT), LSC_spatial_avg, degree = 2)$fit
  
  #LSC_spatial_avg_smooth <- loess.wrapper(c(1:TT), LSC_spatial_avg, span.vals = seq(0.1, 1, by = 0.01))$fit
  #LSC_temporal_avg_smooth <- loess.as(c(1:card_space), LSC_temporal_avg, span.vals = seq(0.1, 1, by = 0.01))$fit
  
  #LSC_spatial_avg_smooth <- loess.as(LSC_spatial_avg ~ c(1:TT), degree = 2)$fit
  
  #field.dim = dim(LSC)
  
  layout(matrix(c(1, 3, 2, 4), byrow = TRUE, ncol = 2), widths = widths, 
         heights = heights)
  par(mar = c(0.5, 5, 2, 0), cex.lab = 2, cex.axis = 2)
  plot(seq_len(card_space), LSC_temporal_avg, main = "", ylab = "", xlab = "", pch = 19, 
       axes = FALSE)
  lines(seq_len(card_space), LSC_temporal_avg_smooth, lwd = 2, col = "red")
  box()
  axis(1, at = pretty(LSC_temporal_avg))
  mtext("bits", 1, line = 3, cex = 1.5)
  mtext("Temporal Average", 3, line = 0.5, cex = 1.5, las = 1)
  
  par(mar = c(5, 5, 0, 0), cex.lab = 2, cex.axis = 2)
  # image(1:ncol(LSC), 1:nrow(LSC), t(LSC), col = col, xlab = 'Time', ylab =
  # 'Space', axes = FALSE)
  image2(LSC, col = col, xlab = "Space", ylab = "Time", axes = FALSE, 
         legend = FALSE)
  box()
  axis(2, at = pretty(seq_len(nrow(LSC)))[-length(pretty(1:nrow(LSC)))], 
       labels = rev(pretty(seq_len(nrow(LSC)))[-1]) )
  axis(1, at = pretty(seq_len(nrow(LSC)))[-length(pretty(1:ncol(LSC)))])
  # mtext('Time', 1, line = 0.5, cex = 1.5) mtext('Space', 2, line = 0.5,
  # cex = 1.5)
  
  par(mar = c(0.5, 0.5, 2, 10))
  make_legend(data = LSC, col = col, side = 4, col.label = lsc.unit)
  
  par(mar = c(5, 0.5, 0, 2), cex.lab = 2, cex.axis = 2)
  plot(LSC_spatial_avg, TT:1,  main = "", ylab = "bits", 
       xlab = "", pch = 19, axes = FALSE)
  lines(LSC_spatial_avg_smooth, TT:1, lwd = 2, col = "red")
  box()
  axis(2, at = pretty(LSC_spatial_avg))
  mtext("Spatial Average", 4, line = 0.5, cex = 1.5)
}

