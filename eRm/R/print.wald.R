`print.wald` <-
function(x,...)
#print method for objects of class "wald" (from waldtest)
{
   #if (!is.null(x$betalab)) {
   #  cat("Warning Message: Item",x$betalab[1],"was not tested due to sum-0 restriction.\n")
   #}
   cat("\nWald test on item level (z-values):\n\n")
   print(round(x$coef.table,3))
   cat("\n")
   invisible(round(x$coef.table,3))
}

