inkei <- function(length=35, squirt=4, thickness=4)
{
   if (thickness <= 3){thickness=4}
   
      cat("          ", rep("_", length), "\n", sep="")
      
      for (t in 1:floor(thickness/2))
      {     cat("          ", rep(" ", length+t), "\\\n", sep="")}
      
      cat("          ", rep(" ", length+floor(thickness/2)), "-o", rep("~", squirt), "\n", sep="")
      
      for (t in 1:(floor(thickness/2)-1))
      {
         cat("          ", rep(" ", length+floor(thickness/2)-t+1), "/\n", sep="")
      }
   
   cat("          ", rep("_", length+1), "/\n", sep="")
   cat("          \\
  .... |...|
  .... |...|
 \\____/____/")
}
#copyright Jim Aloku