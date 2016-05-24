 print.pamrcved <-function(x, ...) {
   cat("Call:\n")
   dput(x$call)
  
   mat <- rbind(threshold = format(round(x$threshold, 3)), nonzero = 
                format(trunc(x$size)), errors = trunc(x$error * nrow(
                                              x$yhat)))
   if(!is.na(x$pvalue.survival[1])){
     mat <- rbind(mat, pvalue=round(x$pvalue.survival,6))
   }
   dimnames(mat) <- list(dimnames(mat)[[1]], paste(1:ncol(mat)))
 
   print(t(mat), quote = FALSE)
   invisible()
 }
