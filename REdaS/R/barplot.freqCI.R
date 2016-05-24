barplot.freqCI <- function(height, percent = TRUE, ...){
  if(class(height) != "freqCI") stop('"height" must be an object of class "freqCI".')
  
  if(percent){
    height$rel_freq <- 100 * height$rel_freq
    height$CIs_low  <- 100 * height$CIs_low 
    height$CIs_high <- 100 * height$CIs_high
  }
  
  defaults <- list(names.arg = height$cat_names,
                   ylim = c(0, max(unlist(height$CIs_high))),
                   col  = gray(.9),
                   main = "Frequencies with Confidence Intervals",
                   xlab = "Categories",
                   ylab = ifelse(percent, "Percent", "Relative Frequencies"))
  these_args <- modifyList(defaults, list(height = as.table(height$rel_freq), lwd=1, ...))
  x_coords_orig <- do.call("barplot", these_args)
  x_coords <- as.numeric(x_coords_orig)
  x_diff <- diff(x_coords)[1L]/4

  y_bars <- cbind(height$CIs_low,height$CIs_high)
  y_ranges <- apply(y_bars, 1, range)
  
  ddd <- modifyList(list(...), list(col="black"))
  
  for(b in seq_len(nrow(y_bars))){
    do.call("segments", modifyList(list(), list(x0=x_coords[b], y0=y_ranges[1L,b], x1=x_coords[b], y1=y_ranges[2L,b], ddd)))
    for(v in seq_len(ncol(y_bars))){
      do.call("segments", modifyList(list(), list(x0=x_coords[b]-x_diff, y0=y_bars[b,v], x1=x_coords[b]+x_diff, y1=y_bars[b,v], ddd)))
    }
  }
  
  invisible(as.numeric(x_coords_orig))
  
}
