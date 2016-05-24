profileCIs <- function(x, y, wt, B, B.inf, side, alpha, l.max, step, useAstar, penalized) {

   p <- ncol(x) # p covariates
   J <- ncol(y) # J categories

   tol1 <- 1e-03
   if (!penalized) tol2 <- 1e-02 else tol2 <- 1e-04

   l.null <- l.max - (0.5 * qchisq(1 - alpha, df = 1))
   endpoints <- matrix(data = 0, nrow = p * J, ncol = 1)

   B.keep <- B
   B.inf.keep <- B.inf

   for (i in 1:(p * J)) {
   
      B <- B.keep
      B.inf <- B.inf.keep
   
      vecB <- vec(t(B))

      if (!penalized) {
        vecB.inf <- vec(t(B.inf))
         if ((side == 1 && vecB.inf[i] == Inf) | (side == -1 && vecB.inf[i] == -Inf)) {
            endpoints[i] <- vecB.inf[i]
            next
         }
      }

      l <- l.max
      iter <- 1 # Only incremented for MLEs; maximum number of iterations is 300
      
      while (l >= l.null & iter <= 300) {
         vecB.last <- vecB
         vecB[i] <- vecB[i] + (step * side)
         B <- matrix(data = vecB, nrow = p, ncol = J, byrow = TRUE)
         if (!penalized) {
            temp <- profile(x, y, wt, B, i, tol2, useAstar, penalized); B <- temp$B; l <- temp$l
         } else {
         temp <- profile(x, y, wt, B, i, tol1, useAstar, penalized); B <- temp$B; l <- temp$l
         }
         vecB <- vec(t(B))
         if (!penalized) iter <- iter + 1
      }
      
      vecB.plus <- vecB.last
      vecB.minus <- vec(t(B))
      
      converge <- FALSE; iter <- 1

      while (!converge && iter < 30) {
         vecB <- 0.5 * (vecB.plus + vecB.minus)
         B <- matrix(data = vecB, nrow = p, ncol = J, byrow = TRUE)
         if (!penalized) {
            temp <- profile(x, y, wt, B, i, tol2, useAstar, penalized); B <- temp$B; l <- temp$l
         } else {
         temp <- profile(x, y, wt, B, i, tol1, useAstar, penalized); B <- temp$B; l <- temp$l
         }
         converge <- (sum(abs(vecB.minus - vecB.plus)) <= tol2 && abs(l - l.null) <= tol2)
         if (l < l.null) vecB.minus <- vec(t(B)) else vecB.plus <- vec(t(B))
         iter <- iter + 1
      }
      
      vecB <- 0.5 * (vecB.minus + vecB.plus); B <- matrix(data = vecB, nrow = p, ncol = J, byrow = TRUE)
      
      if (!penalized) {
         temp <- profile(x, y, wt, B, i, tol2, useAstar, penalized); B <- temp$B; l <- temp$l
      } else {
         temp <- profile(x, y, wt, B, i, tol1, useAstar, penalized); B <- temp$B; l <- temp$l
      }
      
      endpoints[i] <- vecB[i]

      temp <- SASalgo(x, y, wt, B, l.null, tol1, i, side, useAstar, penalized); B <- temp$B ; converge <- temp$converge
      
      if (converge) endpoints[i] <- vecB[i] else is.na <- endpoints[i]

    }

    endpoints <- matrix(data = endpoints, nrow = p, ncol = J, byrow = TRUE)

}
