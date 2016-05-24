corPlotter <- function(x, y, dat, write = FALSE, plot.format = NULL, 
                       yname = NULL){
  if(is.null(yname)){
    yname = y
  }
  p <- ggplot(data = dat, aes_string(x = x, y = y)) +
    geom_point() +
    geom_smooth(method = "lm", colour = "red", alpha = 0.3, lwd = 2) +
    theme_bw() +
    theme(text = element_text(size=25)) +
    labs(y = yname, x = "Mean number of alleles")
  if(write){
    if(plot.format == "png"){
      ggsave(paste(y, "-corPlot.", plot.format, sep = ""), plot = p, 
             width = 8, height = 7, units = "in")
    } else if(plot.format == "eps"){
      ggsave(paste(y, "-corPlot.", plot.format, sep = ""), plot = p, 
             device=cairo_ps, width = 8, height = 7, units = "in")
    } else {
      stop("Plots not written to file. Invalid plot format given.")
    }
  }
  return(p)
}