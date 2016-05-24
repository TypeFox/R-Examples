summary.fa <- function(object, ...){

   nCharObj <- nchar(object)
   cat("Summary of fa object\n")
   cat("---------------\n")
   cat("Sequences      :",length(object),"\n")
   cat("Average length :",mean(nCharObj),"\n")
   cat("Minimum length :",min(nCharObj),"\n")
   cat("Maximum length :",max(nCharObj),"\n")
   invisible(object)

} 
