plot.cv.biglasso <- function(x, log.l = TRUE, type = c("cve", "rsq", "scale", 
                                                       "snr", "pred", "all"), 
                             selected = TRUE, vertical.line = TRUE, col = "red", ...) {
  # inherits cv.ncvreg
  class(x) <- 'cv.ncvreg'
  plot(x = x, log.l = log.l, type = type, selected = selected, 
       vertical.line = vertical.line, col = col, ...)
}
