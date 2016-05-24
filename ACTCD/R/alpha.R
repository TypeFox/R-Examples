alpha <-
function(K){

  M <- 2^K
  alpha <- matrix(0, M, K)
  len <- 1
  for (i in 1:K)
  {
    com <- t(combn(K,i))
    for (j in 1:dim(com)[1])
    {
      alpha[len+j,com[j,]] <- 1
    }
    len <- len+dim(com)[1]
  }
  return(alpha)
}
