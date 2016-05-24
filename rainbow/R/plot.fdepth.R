plot.fdepth <- function(x, show.legend = TRUE, pos.legend = "bottomleft", ...)
{
  plot(x$data, col = gray(.7), ...)
  lines(x$output$median, col = "red", ...)
  lines(x$output$mtrim, col = "blue", ...)
  if (show.legend){
      legend(pos.legend, legend = c("Trimmed mean", "Median"), lty = 1, 
             col = c("blue", "red"), ...) 
  } 
}

