# @param x X coordinates of boxes to be plotted.
# @param y Y coordinates of boxes to be plotted.
# @param color Color of each box.
# @param n_xboxes Width of plot.
# @param n_yboxes Height of plot.
# @param text_ Text in each box. (Optional)
# @param text_col Color of text in each box. (Optional)

#' @importFrom graphics par plot rect text
plotbox <- function(x, y, color, border_col, n_xboxes, n_yboxes,
                    text_ = NULL, text_col, text_cex = 1, text_font = NULL){
  # Store old par
  old_xaxs <- par()$xaxs
  old_yaxs <- par()$yaxs

  # Make X and Y axis have precise ends
  par(xaxs = "i", yaxs = "i")

  plot(0,0, col = "white", xlim = c(0, n_xboxes), ylim = c(0, n_yboxes),
       bty = "n", xaxt = "n", yaxt = "n", ann = FALSE, asp = 1)

  for(i in 1:length(x)){
    rect(x[i], y[i], x[i] + 1, y[i] + 1,
         col = color[i],
         border = border_col[i])
    if(!is.null(text_)){
      text(x[i] + .5, y[i] + .5, text_[i], col = text_col[i], cex = text_cex,
           family = text_font)
    }
  }

  # Restore old par
  par(xaxs = old_xaxs, yaxs = old_yaxs)

}
