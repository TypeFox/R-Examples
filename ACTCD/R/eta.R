eta <-
function (K,J,Q)
{
  M <- 2^K
  A <- alpha(K)
  tmp <- matrix(NA,M,J)
  for (g in 1:M){ #g is latent pattern
    for (j in 1:J){ #j is item
      tmp[g,j] <- ifelse(all(as.logical(A[g,]^Q[j,])),1,0)
    }
  }
  return(tmp)
}
