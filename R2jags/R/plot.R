plot.rjags <- function(x, display.parallel = FALSE,...){
  plot(x$BUGSoutput, display.parallel = display.parallel,...)
}
