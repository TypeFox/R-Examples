print.ICr <- function(x,...)
{
#print method for objects of class "ICr" (from function "IC")

 cat("\nInformation Criteria: \n")
 print(x$ICtable)
 cat("\n")
}