plot.propagate <- function(x, logx = FALSE, ...)
{
  object <- x
  
  ## plot setup
  par(mfrow = c(2, 1))
  par(mar = c(3, 2, 2, 1))   
  
  ## cut off histogram extreme values
  ## and plot histogram
  resSIM <- object$resSIM
  FILTER <- quantile(resSIM, c(0.01, 0.99), na.rm = TRUE) 
  plotDATA <- resSIM[resSIM > FILTER[1] & resSIM < FILTER[2]]
  plotDATA <- plotDATA[!is.na(plotDATA)]
  
  if (logx) plotDATA <- suppressWarnings(log(plotDATA)) 
    
  ## histogram with fitted skew-normal curve
  HIST <- hist(plotDATA, xlab = "", ylab = "", col = "gray", breaks = 100, 
               main = NULL, cex.main = 1, freq = FALSE,  ...)
  DENS <- density(plotDATA)
  lines(DENS, col = 2, lwd = 2)
  title(main = "Histogram of simulation results with density curve (red)\n and 95% confidence interval (blue)\n", cex.main = 0.7)
  abline(v = object$sim[c(5, 6)], col = "darkblue", lwd = 2)
   
  ## boxplot of simulation result
  boxplot(resSIM, horizontal = TRUE, outline = FALSE, col = "gray", 
          main = "Boxplot of simulation results with 95% confidence interval (blue)\n and first (red)/second (orange) -order moments from Taylor expansion", 
          cex.main = 0.7, ...)
  abline(v = object$sim[c(5, 6)], col = "darkblue", lwd = 2)
  abline(v = object$prop[1], col = "darkred", lwd = 2, lty = 1)
  abline(v = c(object$prop[1] - object$prop[3], object$prop[1] + object$prop[3]) , col = "darkred", lwd = 2, lty = 2)
  abline(v = object$prop[2], col = "darkorange1", lwd = 2, lty = 1)
  abline(v = c(object$prop[2] - object$prop[4], object$prop[2] + object$prop[4]) , col = "darkorange1", lwd = 2, lty = 2)
}
  
  