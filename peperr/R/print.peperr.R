print.peperr <- function(x, ...){ 
   cat("ESTIMATION OF PREDICTION ERROR","\n", "Data split in ", 
      length(x$indices$sample.index), " samples","\n", sep="")
   if (is.list(x$selected.complexity)){
      if ((all.equal(x$complexity,x$selected.complexity)!=TRUE)[1]){
         cat("Selected complexity in full data set: ", "\n")
         print(unlist(x$selected.complexity)) 
         cat("\n", "Selected sample complexity: ", "\n")
         for (i in 1:(length(x$sample.complexity)/2)){
            print(x$sample.complexity[c((i*2)-1, i*2)])
         }
      } else {
          cat("Passed complexity: ", "\n")
          print(x$complexity)
       }
   } else {
      if ((all.equal(x$complexity,x$selected.complexity)!=TRUE)[1]){
         cat("Selected complexity in full data set: ", x$selected.complexity, "\n") 
         cat("Selected sample complexity: ", x$sample.complexity, "\n")
       } else {
          cat("Passed complexity: ", x$sample.complexity, "\n")
       }
   }
}
