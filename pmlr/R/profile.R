profile <- function(x, y, wt, B, i, tol, useAstar, penalized) {

   p <- ncol(x)
   J <- ncol(y)

   ys <- 1 - apply(y, 1, sum)

   e <- matrix(data = 0, nrow = p * J, ncol = 1); e[i,1] <- 1

   vecB <- vec(t(B))

   f1 <- texp(x %*% B); f2 <- f1/(1 + apply(f1, 1, sum)); fref0 <- 1/(1 + apply(f1, 1, sum))

   score <- vec(t((t(x * wt)) %*% (y - f2)))

   if (penalized) {
      if (useAstar) {
         temp <- getV(x, wt, B); A <- temp$A; W <- temp$W; V <- temp$V
      } else {
      temp <- getW(x, wt, B); A <- temp$A; W <- temp$W
      }
   } else {
      A <- getA(x, wt, B)
   }

   l <- apply(t(tln(cbind(f2, fref0)) * cbind(y, ys)) %*% wt, 2, sum)
   if (penalized) l <- l + 0.5 * tln(det(A))

   iter1 <- iter2 <- 1; converge <- FALSE # 'iter1' is used for ONLY for MLEs and 'iter2' used ONLY for PLEs

   while (!converge && iter1 < 6 && iter2 < 4) {
   
      if (penalized) {
         if (det(A) < 1e-05) diag(A) <- diag(A) + 1e-06
         Ainv <- solve(A)
         score <- score + 0.5 * (W %*% vec(Ainv))
         if (det(A) < 1e-05) diag(A) <- diag(A) + 1e-06
         Ainv <- solve(A)
      } else {
         if (det(A) < 1e-02) diag(A) <- diag(A) + 1e-02
         Ainv <- solve(A)
      }

      lambda <- (-(t(e) %*% Ainv %*% score))/Ainv[i,i]
      delta <- Ainv %*% (score + e %*% lambda); delta <- matrix(data = delta, nrow = p, ncol = J, byrow = TRUE)

      iter3 <- 1; increase <- 0

      while (increase == 0 & iter3 < 8) {
         B.temp <- B + delta
         f1 <- texp(x %*% B.temp); f2 <- f1/(1 + apply(f1, 1, sum)); fref0 <- 1/(1 + apply(f1, 1, sum))
         temp <- getW(x, wt, B.temp); A <- temp$A ; W <- temp$W
         l.temp <- apply(t(tln(cbind(f2, fref0)) * cbind(y, ys)) %*% wt, 2, sum); if (penalized) l.temp <- l.temp + 0.5 * tln(det(A))
         increase <- (l.temp > l)
         if (!increase) delta <- delta/2
         iter3 <- iter3 + 1
      }

      B <- B + delta
      f1 <- texp(x %*% B); f2 <- f1/(1 + apply(f1, 1, sum)); fref0 <- 1/(1 + apply(f1, 1, sum))
      score <- vec(t((t(x * wt)) %*% (y - f2)))
      if (useAstar) {
         temp <- getV(x, wt, B); A <- temp$A; W <- temp$W; V <- temp$W
      } else {
      temp <- getW(x, wt, B); A <- temp$A ; W <- temp$W
      }
      l <- apply(t(tln(cbind(f2, fref0)) * cbind(y, ys)) %*% wt, 2, sum); if (penalized) l <- l + 0.5 * tln(det(A))
# Why is the convergence criteria commented out in the penalized case?
      if (!penalized) converge <- any(abs(delta) <= tol)
      if (!converge) {
         if (!penalized) iter1 <- iter1 + 1 else iter2 <- iter2 + 1
      }
   }

   list(B = B, l = l)

}
