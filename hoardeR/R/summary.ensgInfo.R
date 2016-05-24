summary.ensgInfo <- function(object, ...){
  
   objTable <- table(object$GeneType)

   cat("Summary of the ensgInfo object\n")
   cat("------------------------------\n")
   cat("Gene types:\n")
   for(i in 1:length(objTable)){
     cat(paste(names(objTable)[i],":",objTable[i],"\n"))
   }
   cat("------------------------------\n")
   invisible(object)
} 
