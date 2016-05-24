"BEhypergeo" <- function(p, q, N, D, lambda, eps=1E-12, maxit=500) {
   if(lambda < 0) {
     warning("Invalid lambda argument")
     return()
   }
   if(length(N) != p) n <- rep(N[1], p)
   if(length(D) != p) d <- rep(D[1], p)
   k <- -1; err <- Inf; Fpq <- 0; PRIOR <- Inf
   while(1) {
      k <- k + 1
      NK <- sum(sapply(1:p, function(i) lgamma(n[i] + k) + lgamma(n[i])))
      DK <- sum(sapply(1:q, function(i) lgamma(d[i] + k) + lgamma(d[i])))
      TMP <- exp(log(lambda^k) - lgamma(k+1) + NK - DK)
      if(     is.nan(TMP)) break
      if(! is.finite(TMP)) break
      #message("TMP: ",TMP)
      Fpq <- Fpq + TMP
      err <- ((Fpq - PRIOR)/Fpq)^2
      if(err < eps) break;
      if(k == maxit) {
         warning("Maximum iterations reached, result might be unreliable")
         break;
      }
      PRIOR <- Fpq
   }
   return(list(value=Fpq, its=(k+1), error=err))
}
