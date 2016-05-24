numericGradient <- function(f, t0, eps=1e-6, ...) {
   ### numeric gradient of a vector-valued function
   ### f    : function, return Nval x 1 vector of values
   ### t0   : NPar x 1 vector of parameters
   ### return:
   ### NvalxNPar matrix, gradient
   NPar <- length(t0)
   NVal <- length(f0 <- f(t0, ...))
   grad <- matrix(NA, NVal, NPar)
   row.names(grad) <- names(f0)
   colnames(grad) <- names(t0)
   for(i in 1:NPar) {
      t2 <- t1 <- t0
      t1[i] <- t0[i] - eps/2
      t2[i] <- t0[i] + eps/2
      grad[,i] <- (f(t2, ...) - f(t1, ...))/eps
   }
   return(grad)
}