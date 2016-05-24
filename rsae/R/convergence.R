convergence <-
function(object){
   if (!inherits(object, "fitsaemodel")) stop("method convergence is not applicable\n")
   saemodel <- attr(object, "saemodel")
   n <- saemodel$n
   g <- saemodel$g
   p <- saemodel$p
   fixeffnames <- colnames(saemodel$X)
   # retrieve the estimating method
   method <- attr(object, "method") 
   cat("CONVERGENCE REPORT \n")
   # check whether the model converged
   converged <- object$converged
   if (converged != 1){
      cat("NOTE: ALGORITHM DID NOT CONVERGE!\n")
   }
   #----------------------
   # niter and acc specification
   optim <- attr(object, "optim")
   acc <- optim$acc
   niter <- c(optim$niter, 100) # 100 is the max value defined in estimation of "d"
   together <- cbind(niter, acc)
   colnames(together) <- c("niter", "acc")
   rownames(together) <- c("overall loop", "fixeff", "residual var", "area raneff var")
   cat("---\n")
   cat("User specified number of iterations (niter) and \nnumeric precision (acc):\n")
   cat("\n")
   print.default(format(together, digits = 1), print.gap = 2, quote = FALSE)
   #----------------------
   # used iters
   iters <- optim$usediter
   if (dim(iters)[1] == 1){
      iters <- as.matrix(iters)
   }else{
      #remove the zero entires in iters
      iters <- iters[rowSums(iters) != 0 ,]
   }
   colnames(iters) <- c("fixeff", "residual var", "area raneff var")
   rownames(iters) <- paste(seq(1, dim(iters)[1]))
   cat("---\n")
   cat(paste("Number of runned EE-specific iterations in each \ncall (given the user-defined specs), reported for \neach of the ",dim(iters)[1], "overall iterations separately: \n"))
   cat("\n")
   print.default(format(iters, digits = 1), print.gap = 2, quote = FALSE)
}

