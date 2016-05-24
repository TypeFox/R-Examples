
crosshair <- function(x, ...) UseMethod("crosshair")

crosshair.default <- function(x, correction = 0.5, level = 0.95, method = "wilson", 
                      xlim = c(0,1), ylim = c(0,1), length = 0.1, pch = 1, 
                      add = FALSE, suppress = TRUE, ...){
if(suppress){x <- suppressWarnings(  x <- madad(x, correction = correction, level = level, 
             method = method))
             }else{
               x <-   x <- madad(x, correction = correction, level = level, 
             method = method)
             }
  if(!add){plot(x$fpr$fpr, x$sens$sens, xlim = xlim, ylim = ylim, pch = pch, 
                xlab = "False Positive Rate", ylab = "Sensitivity", ...)}
  if(add){points(x$fpr$fpr, x$sens$sens, pch = pch, ...)}
  arrows(x$fpr$fpr.ci[,1], x$sens$sens,  x$fpr$fpr.ci[,2], x$sens$sens, angle = 90, code = 3, length = length, ...)
  arrows(x$fpr$fpr, x$sens$sens.ci[,1],  x$fpr$fpr, x$sens$sens.ci[,2],  angle = 90, code = 3, length = length, ...)
  return(invisible(NULL))
  }
