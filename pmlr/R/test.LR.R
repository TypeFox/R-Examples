test.LR <- function(x, y, wt, B, h0, penalized) {

   p <- nrow(B) # p covariates
   J <- ncol(B) # J categories

   tol <- 1e-6
   ys <- 1 - apply(y,1,sum)

# under H_a
# is this the likelihood from fit$lstar.max?

   f1 <- texp(x %*% B); f2 <- f1/(1 + apply(f1, 1, sum)); fref0 <- 1/(1 + apply(f1, 1, sum))
   temp <- getW(x, wt, B); A <- temp$A; W <- temp$W
   la <- apply(t(tln(cbind(f2, fref0)) * cbind(y, ys)) %*% wt, 2, sum)
   if (penalized) la <- la + 0.5 * tln(det(A))

# calculate the likelihhod under H_0

   if (h0 == 1) {ntests <- p * J;   size <- J;  L <- 1}
   if (h0 == 2) {ntests <- p;       size <- 1;  L <- J}
   if (h0 == 3) {ntests <- p;       size <- 1;  L <- J - 1}

   l0 <- pvalue <- matrix(data=0, nrow=ntests, ncol=1)

   for (i in 1:ntests) {
      iter1 <- 0
      if (h0 == 1) {e <- matrix(data=0, nrow=p * J, ncol=L); e[i,1] <- 1}
      if (h0 == 2) {C1 <- matrix(data=0, nrow=L, ncol=p * J); C1[1:L, (((i - 1) * J) + 1):(i * J)] <- diag(L); e <- t(C1)}
      if (h0 == 3) {
         C1 <- matrix(data=0, nrow=L, ncol=p * J)
         C1[1:L, (((i - 1) * J) + 1):((i * J) - 1)] <- diag(L)
         for (m in 1:(J - 1)) {C1[m, ((i - 1) * J) + 1 + m] <- (-1)}
         e <- t(C1)
      }

      B0 <- matrix(data = 0, nrow = p, ncol = J, byrow = TRUE)
      converge <- FALSE

      while (!converge & iter1 <= 20) {
         f1 <- texp(x %*% B0); f2 <- f1/(1 + apply(f1, 1, sum)); fref0 <- 1/(1 + apply(f1, 1, sum))
         temp <- getW(x, wt, B0); A <- temp$A; W <- temp$W
         l0[i,1] <- apply(t(tln(cbind(f2, fref0)) * cbind(y, ys)) %*% wt, 2, sum)
         if (penalized) l0[i,1] <- l0[i,1] + 0.5 * tln(det(A))

# why is this 1e-02 here?
         if (det(A) < 0.01) diag(A) <- diag(A) + 0.01
         Ainv <- solve(A)

         score <- vec(t((t(x * wt)) %*% (y - f2)))
         if (penalized) score <- score + 0.5 * (W %*% vec(Ainv))

         lambda <- -solve(t(e) %*% Ainv %*% e) %*% (t(e) %*% Ainv %*% score)
         delta <- Ainv %*% (score + e %*% lambda); delta <- matrix(data = delta, nrow = p, ncol = J, byrow = TRUE)

         iter2 <- 1; increase <- FALSE;
         while (!increase & iter2 < 5) {
            B0.temp <- B0 + delta
            f1 <- texp(x %*% B0.temp); f2 <- f1/(1 + apply(f1, 1, sum)); fref0 <- 1/(1 + apply(f1, 1, sum))
            temp <- getW(x, wt, B0.temp); A <- temp$A; W <- temp$W
            l.temp <- apply(t(tln(cbind(f2, fref0)) * cbind(y, ys)) %*% wt, 2, sum)
            if (penalized) l.temp <- l.temp + 0.5 * tln(det(A))
            increase <- (l.temp > l0[i,1])
            if (!increase) delta <- delta/2
            iter2 <- iter2 + 1
         }

         B0 <- B0 + delta
         f1 <- texp(x %*% B0); f2 <- f1/(1 + apply(f1, 1, sum)); fref0 <- 1/(1 + apply(f1, 1, sum))
         temp <- getW(x, wt, B0); A <- temp$A; W <- temp$W
         l0[i,1] <- apply(t(tln(cbind(f2, fref0)) * cbind(y, ys)) %*% wt, 2, sum)
         if (penalized) l0[i,1] <- l0[i,1] + 0.5 * tln(det(A))
         converge <- max(abs(delta)) < tol
         iter1 <- iter1 + 1
      }

   }

   statistic <- matrix(data = LR <- 2 * (la - l0), nrow = p, ncol = size, byrow = TRUE)
   pvalue <- pchisq(statistic, L, lower.tail = FALSE)

   list(statistic = statistic, pvalue = pvalue)

}
