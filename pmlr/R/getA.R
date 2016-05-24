# A = pJ x pJ Fisher information matrix (based on the 2nd derivatives of the log-likelihood with respect to B)

getA <- function(x, wt, B) {

   p <- nrow(B) # p covariates
   J <- ncol(B) # J categories

   f1 <- texp(x %*% B)
   f2 <- f1/(1 + apply(f1, 1, sum))

   A <- matrix(data = 0, nrow = p*J, ncol = p*J)

   ij <-1
   while (ij <= J) {
      ik <- 1
      while (ik <= J) {
         if (ij == ik) {
            s <- f2[,ij] * (1 - f2[,ik]); cc1 <- 1
         } else {
            s <- f2[,ij] * f2[,ik]; cc1 <- -1
         }
         sx <- sqrt(s) * x * sqrt(wt)
         A <- A + cc1 * (crossprod(sx) %x% (cbind(diag(J)[,ij]) %*% rbind(diag(J)[ik,])))
         ik <- ik + 1
      }
      ij <- ij + 1
   }

   A

}
