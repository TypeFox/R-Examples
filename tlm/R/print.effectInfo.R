print.effectInfo <-
function(x, ...)
 {  
  if (!inherits(x, "effectInfo"))
    stop("argument 'x' must be of class 'effectInfo'")
  if (attr(x, "modeltype") > 0)
   {
   	cat("\n")
   	cat("The effect of X on Y can be summarized with a single number as follows:\n\n")
   	cat(" - Change in X:", x$Xincrease, "\n")
   	cat(" - Type of effect on Y:", x$effecttype, "\n")
   	cat(" - Effect size:", x$effectsize, "\n\n")
   	cat("   beta coefficient estimate:\n")
   	print(x$beta, ...)
   	cat("\nFurther details can be obtained using effect().\n\n")
    } else cat("\nThe effect of X on Y cannot be summarized with a single number.\nIts behavior can be explored using effect().\n\n")
 }
