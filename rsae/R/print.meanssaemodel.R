print.meanssaemodel <-
function(x, digits=4, ...){
   cat("Robustly Estimated/Predicted Area-Level Means:\n")
   hasmspe <- attr(x, "mspe")
   if (is.null(hasmspe)){
      all <- cbind(x$raneff, x$fixeff, x$means)
      colnames(all) <- c("raneff", "fixeff", "predicted mean")
   }else{
      all <- cbind(x$raneff, x$fixeff, x$means, x$mspe)
      colnames(all) <- c("raneff", "fixeff", "area mean", "MSPE")
     }
   print.default(format(all, digits = digits), print.gap = 2, quote = FALSE)
   if (!is.null(hasmspe)){
   cat(paste("(MSPE: ", hasmspe," boostrap replicates)\n", sep=""))
   }
}

