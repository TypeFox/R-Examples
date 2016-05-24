SASalgo <- function(x, y, wt, B, l.null, tol, i, side, useAstar = NULL, penalized) {

   p <- nrow(B) # p covariates
   J <- ncol(B) # J categories

   vecB <- vec(t(B))

   e <- matrix(data = 0, nrow = p * J, ncol = 1); e[i,1] <- 1

   converge <- FALSE; iter1 <- 1

   while (!converge & iter1 < 25) {

      f1 <- texp(x %*% B); f2 <- f1/(1 + apply(f1, 1, sum)); fref0 <- 1/(1 + apply(f1, 1, sum)); ys <- 1 - apply(y, 1, sum)
      score <- vec(t((t(x * wt)) %*% (y - f2)))
      l <- apply(t(tln(cbind(f2, fref0)) * cbind(y, ys)) %*% wt, 2, sum)

      if (!penalized) {
         A <- getA(x, wt, B)
         if (det(A) < 1e-05) diag(A) <- diag(A) + 1e-05
         Ainv <- solve(A)
      } else {
         if (useAstar) {
            temp <- getV(x, wt, B); A <- temp$A; W <- temp$W; V <- temp$V
         } else {
         temp <- getW(x, wt, B); A <- temp$A; W <- temp$W
         }
         l <- l + 0.5 * tln(det(A))
# Why is this adjustment to the Fisher information matrix made AFTER penalizing the likelihood instead of before?
         if (det(A) < 1e-05) diag(A) <- diag(A) + 1e-06
         Ainv <- solve(A)
         score <- score + 0.5 * (W %*% vec(Ainv))
         if (useAstar) {
            A <- getAstar(p, J, A, W, V)
            if (det(A) < 1e-05) diag(A) <- diag(A) + 1e-06
            Ainv <- solve(A)
         }
      }

      if (l < l.null + 0.5 * (t(score) %*% Ainv %*% score)) {
         lambda <- (-(t(e) %*% Ainv %*% score))/Ainv[i,i]
         delta <- Ainv %*% (score + e %*% lambda)
      } else {
         lambda <- side * sqrt(2 * (l - l.null + 0.5 * t(score) %*% Ainv %*% score)/Ainv[i,i])
         delta <- Ainv %*% (score + e %*% lambda)

         if (!penalized) {
            inside <- FALSE; iter2 <- 1
            while (!inside & iter2 <= 11) {
               vecB.temp <- vecB + delta
               B.temp <- matrix(data = vecB.temp, nrow = p, ncol = J, byrow = TRUE)
               f1 <- texp(x %*% B.temp); f2 <- f1/(1 + apply(f1, 1, sum)); fref0 <- 1/(1 + apply(f1, 1, sum))
               l <- apply(t(tln(cbind(f2, fref0)) * cbind(y,ys)) %*% wt, 2, sum)
               inside <- l > l.null
               if (!inside) {delta <- delta/2; iter2 <- iter2 + 1}
            }
         }
      }
      
      vecB <- vecB + delta; B <- matrix(data = vecB, nrow = p, ncol = J, byrow = TRUE)

      diff <- l - l.null
      if (!penalized) normdelta <- max(abs(delta)) else normdelta <- t(score + e %*% lambda) %*% Ainv %*% (score + e %*% lambda)
      converge <- (abs(diff) <= tol & normdelta <= tol)
      if (!converge) iter1 <- iter1 + 1

    }

# Should we be using this?
    if (!penalized) converge <- abs(diff) < 0.01
    return(list(B = B, converge = converge))

}
