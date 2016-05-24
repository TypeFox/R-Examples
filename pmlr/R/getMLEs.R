getMLEs <- function(x, y, wt) {

   p <- ncol(x) # p covariates
   J <- ncol(y) # J categories

   B <- matrix(data = 0, nrow = p, ncol = J)
   vecB <- vec(t(B))
   
   Ainv <- diag(p * J)
   ratio <- rep(1, length = p * J)

   for (iter in 1:49) {
      f1 <- texp(x %*% B); f2 <- f1/(1 + apply(f1, 1, sum)); fref0 <- 1/(1 + apply(f1, 1, sum))
      score <- vec(t((t(x * wt)) %*% (y - f2)))
      A <- getA(x, wt, B)
      Ainv.old <- Ainv
      if (inherits(try(solve(A), silent = TRUE), "try-error")) break else Ainv <- solve(A)
      Ainv.new <- Ainv
      step <- Ainv %*% score
      vecB <- vecB + step
      B <- matrix(data = vecB, nrow = p, ncol = J, byrow = TRUE)
      ratio <- diag(Ainv.new)/diag(Ainv.old)
   }

   inf <- diag(Ainv) > 500 | ratio > 2

   B <- matrix(data = 0, nrow = p, ncol = J)
   vecB <- vec(t(B))
   
   ys <- 1 - apply(y, 1, sum)

   iter1 <- 1; stop <- FALSE; tol <- 1e-05

   while (!stop) {
   
      f1 <- texp(x %*% B); f2 <- f1/(1 + apply(f1, 1, sum)); fref0 <- 1/(1 + apply(f1, 1, sum))
      score <- vec(t((t(x * wt)) %*% (y - f2)))
      l <- apply(t(tln(cbind(f2, fref0)) * cbind(y, ys)) %*% wt, 2, sum)
      A <- getA(x, wt, B)
      if (inherits(try(solve(A), silent = TRUE), "try-error")) break else Ainv <- solve(A)

      step <- Ainv %*% score
      iter2 <- 1; increase <- FALSE

      while (!increase && iter2 <= 10) {
         vecB.temp <- vecB + step
         B.temp <- matrix(data = vecB.temp, nrow = p, ncol = J,byrow=TRUE)
         f1 <- texp(x %*% B.temp); f2 <- f1/(1 + apply(f1, 1, sum)); fref0 <- 1/(1 + apply(f1, 1, sum))
         l.temp <- apply(t(tln(cbind(f2, fref0)) * cbind(y,ys)) %*% wt, 2, sum)
         increase <- l.temp > l
         if (!increase) { step <- step/2; iter2 <- iter2 + 1 }
      }

      vecB <- vecB + step; B <- matrix(data = vecB, nrow = p, ncol = J, byrow = TRUE)
      vecB.inf <- ifelse(inf, vecB/0, vecB); B.inf <- matrix(data = vecB.inf, nrow = p, ncol = J, byrow = TRUE)
#     stop <- max(abs(step)[!inf]) <= tol || iter1 == 15
      stop <- iter1 == 15
      if (!stop) iter1 <- iter1 + 1

   }

   f1 <- texp(x %*% B); f2 <- f1/(1 + apply(f1, 1, sum)); fref0 <- 1/(1 + apply(f1, 1, sum))
   A <- getA(x, wt, B)
   l.max <- apply(t(tln(cbind(f2, fref0)) * cbind(y, ys)) %*% wt, 2, sum)
   if (inherits(try(solve(A), silent = TRUE), "try-error")) break else Ainv <- solve(A)
   
   fit <- list(B = B, B.inf = B.inf, Ainv = Ainv, l.max = l.max, iter1 = iter1)

   fit

}
