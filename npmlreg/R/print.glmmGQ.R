"print.glmmGQ" <-
function(x,digits=max(3,getOption('digits')-3), ...)
{
  cat("\nCall: ",deparse(x$call),"\n\n")
  cat("Coefficients:\n")
  print.default(format(x$coefficients, digits = digits), print.gap = 2, quote = FALSE)
  if (x$family$family=='gaussian'&& x$Misc$lambda<=1/length(x$masses)){
      cat("\nComponent distribution - MLE of sigma:\t  ",format(signif(x$sdev$sdev,digits)),"\n")# 03/09/07
  }
  if (x$family$family=='Gamma'&& x$Misc$lambda<=1/length(x$masses)){ #print shape only if it is constant over components
       cat("\nComponent distribution - MLE of shape parameter:\t  ",format(signif(x$shape$shape,digits)),"\n")# 03/09/07
  }
  cat("Random effect distribution - standard deviation:\t  ", format(x$rsdev),"\n\n") # 03/09/07
  cat("-2 log L:\t   ",format(round(x$disparity,digits=1)),"\n")
  invisible(x)
}

