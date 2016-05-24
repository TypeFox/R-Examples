
print.dBox <- function(x, digits = max(3, getOption("digits") - 3), printCall = TRUE, ...)
{
   cat("Box-like desirability function\n")
   cat("\nCall: ", deparse(x$call), "\n\n", sep = "")
   
   cat("Non-informative value:", round(x$missing, digits), "\n")
   if(!is.null(x$tol)) cat("tolerance:", round(x$tol, digits), "\n")

   invisible(x)
}

print.dMax <- function(x, digits = max(3, getOption("digits") - 3), printCall = TRUE, ...)
{
   cat("Larger-is-better desirability function\n")
   cat("\nCall: ", deparse(x$call), "\n\n", sep = "")
   
   cat("Non-informative value:", round(x$missing, digits), "\n")
   if(!is.null(x$tol)) cat("tolerance:", round(x$tol, digits), "\n")

   invisible(x)
}

print.dMin <- function(x, digits = max(3, getOption("digits") - 3), printCall = TRUE, ...)
{
   cat("Smaller-is-better desirability function\n")
   cat("\nCall: ", deparse(x$call), "\n\n", sep = "")
   
   cat("Non-informative value:", round(x$missing, digits), "\n")
   if(!is.null(x$tol)) cat("tolerance:", round(x$tol, digits), "\n")
   
   invisible(x)
}

print.dTarget <- function(x, digits = max(3, getOption("digits") - 3), printCall = TRUE, ...)
{
   cat("Target-is-best desirability function\n")
   cat("\nCall: ", deparse(x$call), "\n\n", sep = "")
   
   cat("Non-informative value:", round(x$missing, digits), "\n")

   if(!is.null(x$tol)) cat("tolerance:", round(x$tol, digits), "\n")
   invisible(x)
}

print.dArb <- function(x, digits = max(3, getOption("digits") - 3), printCall = TRUE, ...)
{
   cat("Arbitrary desirability function\n")
   cat("\nCall: ", deparse(x$call), "\n\n", sep = "")
   
   cat("Non-informative value:", round(x$missing, digits), "\n")
   if(!is.null(x$tol)) cat("tolerance:", round(x$tol, digits), "\n")

   invisible(x)
}


print.dCategorical <- function(x, digits = max(3, getOption("digits") - 3), printCall = TRUE, ...)
{
   cat("Desirability function for categorical data\n")
   cat("\nCall: ", deparse(x$call), "\n\n", sep = "")
   
   cat("Non-informative value:", round(x$missing, digits), "\n")
   if(!is.null(x$tol)) cat("tolerance:", round(x$tol, digits), "\n")

   invisible(x)
}


print.dOverall <- function(x, digits = max(3, getOption("digits") - 3), printCall = TRUE, ...)
{
   cat("Combined desirability function\n")
   cat("\nCall: ", deparse(x$call), "\n\n", sep = "")

   for(i in seq(along = x$d))
   {
      cat("----\n")
      print(x$d[[i]], printCall = FALSE)
   }     
   invisible(x)
}
