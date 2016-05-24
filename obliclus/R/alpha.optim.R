alpha.optim <-
function(A, T, G, cluster, info, maxit=1000){
  alpha <- info$alpha

  for(i in 1:maxit){
    X <- T - alpha * G
    T.proj <- X %*% sqrt(solve(diag(diag(t(X) %*% X))))
    h.alpha <- criteria(A, T.proj, cluster, info)
    h.0 <- criteria(A, T, cluster, info)
    if(h.alpha < h.0)
      break
    else
      alpha <- alpha / 2
    if(i == maxit){
      alpha <- 0
    }
  }

  return(alpha)
}
