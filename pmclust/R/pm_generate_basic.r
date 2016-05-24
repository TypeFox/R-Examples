### This file contains a simple data generation.

### N.K.spmd: integer[K], number of elements of each cluster.
generate.basic.spmd <- function(N.allspmds, N.spmd, N.K.spmd, N, p, K){
  data.simu <- NULL
  data.class <- NULL
  data.n.class <- rep(0, K)
  for(i.k in 1:K){
    tmp.n.k <- N.K.spmd[i.k]

    if(tmp.n.k > 0){
      tmp.data.simu <- NULL
      for(i.p in 1:p){
        mean <- i.k * 2 + i.p + 3
        sd <- sqrt(1 / mean)
        tmp.data.simu <- cbind(tmp.data.simu, rnorm(tmp.n.k, mean, sd))
      }
      data.simu <- rbind(data.simu, tmp.data.simu)
      data.class <- c(data.class, rep(i.k, tmp.n.k))
      data.n.class[i.k] <- tmp.n.k
    }
  }

  ret <- list(K = K, p = p, N = N, N.allspmds = N.allspmds,
              N.spmd = N.spmd, N.K.spmd = N.K.spmd,
              X.spmd = data.simu, CLASS.spmd = data.class,
              N.CLASS.spmd = data.n.class)
  ret
} # End of generate.basic.spmd().

generate.basic <- generate.basic.spmd
