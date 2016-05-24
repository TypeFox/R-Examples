"print.summary.bounds" <-
function(x, digit = 5, ...)
{
   z <- x
   if (!inherits(z, "summary.bounds")) 
     stop("'x' must inherit from class \"summary.bounds\"")
   rownames(z$bounds) <- rownames(z$bounds, do.NULL = FALSE, prefix = "")
   cat("\nLan-DeMets bounds for a given spending function \n", "\nn = ", z$n, "\nOverall alpha: ", z$oalpha, "\n")
   if (z$type==1){
      cat("\nType: One-Sided Bounds", "\nalpha: ", z$alpha, "\nSpending function:", z$spending, "\n", "\nBoundaries:\n")
   if (ncol(z$bounds)==5) 
      print.default(z$bounds[,-2], digits = 5, quote = FALSE, print.gap = 2)
   else
      print.default(z$bounds[,-3], digits = 5, quote = FALSE, print.gap = 2)
   }
   else{
      if (z$type==2){
 cat("\nType: Two-Sided Symmetric Bounds", "\nLower alpha: ", z$alpha[1], "\nUpper alpha: ", z$alpha[1], "\nSpending function: ", z$spending, "\n")
      }                                               
      if (z$type==3){
 cat("\nType: Two-Sided Asymmetric Bounds", "\nLower alpha: ", z$alpha[1], "\nSpending function for the lower boundary: ", z$spending[1], "\nUpper alpha: ", z$alpha[2], "\nSpending function for the upper boundary: ", z$spending[2], "\n")
      }
      cat("\nBoundaries:\n")
      print.default(z$bounds, digits = digit, quote = FALSE, print.gap = 2)
   }
}

