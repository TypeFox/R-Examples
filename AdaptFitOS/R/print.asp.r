########## R-function: print.spm ##########


print.asp <- function(x,...)
{
   object <- x

   cat("\n This is an asp2() fit object; a list with components named:\n\n")

   cat("     ",names(object),"\n\n")

  summary(x)
}










