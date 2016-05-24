"print.summary.drift" <-
function(x, digit = 5, ...)
{
   z <- x
   if (!inherits(z, "summary.drift")) 
     stop("'x' must inherit from class \"summary.drift\"")
   rownames(z$bounds) <- rownames(z$bounds, do.NULL = FALSE, prefix = "")
   cat("\nLan-DeMets method for group sequential boundaries \n", "\nn = ", z$n, "\n")
   cat("\nBoundaries: \n") 
   print.default(z$bounds, digits = digit, quote = FALSE, print.gap = 2)
   if ((z$type==1)|(z$type==2)){
      if (z$type==1){
        cat("\nPower : ", z$power, "\n","\nDrift: ", z$drift, "\n")
      }
      if (z$type==2){
        cat("\nDrift parameters: ", z$drift, "\n")
      }
      rownames(z$bounds1) <- rownames(z$bounds1, do.NULL = FALSE, prefix = "")
      print.default(z$bounds1, digits = digit, quote = FALSE, print.gap = 2)
   }
   if (z$type==3){
     low <- z$interval$lower.limit
     up <- z$interval$upper.limit
     cat("\nConfidence interval at the end of the trial: \n", "\nConfidence level: ", z$level, "\nLast Z value: ", z$fzvalue, "\n", 100*z$level, "% confidence interval: (", low, ",", up, ")\n") 
   }
}

