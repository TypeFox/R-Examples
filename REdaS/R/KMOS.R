KMOS <- function(x, use = c("everything", "all.obs", "complete.obs", "na.or.complete", "pairwise.complete.obs")){
  
  this_call <- match.call()
  use <- match.arg(use)

  cormat  <- cor(x = x, use = use)
  pcormat <- solve(cormat)

  n <- nrow(x)
  k <- ncol(x)
  seq_k <- seq_len(k)

  for(j in seq_k){
    for(i in seq_k){
      if(i == j){
        next
      } else {
        pcormat[i,j] <- - pcormat[i, j] / sqrt(pcormat[i, i] * pcormat[j, j])
      }
    }
  }

  MSA <- unlist(lapply(seq_k, function(i){
    sum((cormat[i,-i])^2) / ( sum((cormat[i,-i])^2) + sum((pcormat[i,-i])^2) )
  }))
  names(MSA) <- colnames(x)
  
  KMO <- sum((cormat[!diag(k)])^2) / ( sum((cormat[!diag(k)])^2) + sum((pcormat[!diag(k)])^2) )
  
  res <- structure(
    list("call"    = this_call,
         "cormat"  = cormat,
         "pcormat" = pcormat,
         "n"       = n,
         "k"       = k,
         "MSA"     = MSA,
         "KMO"     = KMO),
    "class" = "MSA_KMO"
  )
  return(res)

}
