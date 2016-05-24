summary.xmlImport <- function(object, ...){

   rowCount <- sapply(object,nrow)
   cat("Summary of xmlImport object\n")
   cat("---------------------------\n")
   cat("Sequences    :",length(object),"\n")
   cat("Min hits     :",min(rowCount),"\n")
   cat("Average hits :",mean(rowCount),"\n")
   cat("Median hits  :",median(rowCount),"\n")
   cat("Max hits     :",max(rowCount),"\n")
   invisible(object)

} 
