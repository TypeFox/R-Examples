summary.saemodel <-
function(object, ...){
   cat("Model summary: \n")
   print(attr(object, "call"))
   cat("---\n")
   nsize <- object$nsize
   cat("No. of areas: ", length(nsize), "\n")
   cat("No. of obs.: ", sum(nsize), "\n")
   if (length(unique(nsize)) == 1){
      cat("Balanced data, each area has ", nsize[1], " units \n")
   } 
   else{
      cat("Smallest area: ", min(nsize), " units \n")
      cat("Largest area: ", max(nsize), " units \n")
   }
}

