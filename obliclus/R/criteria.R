criteria <-
function(A, T, cluster, info){
  w <- info$oblimin.index
  N.var <- dim(A)[1]
  N.fac <- dim(A)[2]
  N.cluster <- info$N.cluster
  L <- A %*% solve(t(T))

  cr.l <- cr.r <- 0

  ##Value of rotation criteria
  if (info$method == "oblimin") {
    Q <- matrix(0, N.fac, N.fac)
    for (i in 1:N.fac) {
      for (j in 1:N.fac) {
        Q[i,j] <- sum(L[,i]^2 * L[,j]^2) - w / N.var * sum(L[,i]^2) * sum(L[,j]^2)
      }
    }
    diag(Q) <- 0
    cr.l <- sum(Q)
  }
  else if (info$method == "geomin") {
    k <- ncol(A)
    p <- nrow(A)
    delta <- info$geomin.par
    L2 <- L^2 + delta
    pro <- exp(rowSums(log(L2))/k)
    cr.l <- sum(pro)
  }

  ##Value of k-means criterion
  L2.mean <- matrix(0, N.cluster, N.fac)
  for(j in 1:N.var)
    L2.mean[cluster[j],] <- L2.mean[cluster[j],] + L[j,]^2

  N.clust.in <- numeric(N.cluster)
  for (k in 1:N.cluster)
    N.clust.in[k] <- length(which(cluster == k))
  L2.mean <- L2.mean / N.clust.in

  for(j in 1:N.var)
    cr.r <- cr.r + sum((L[j,]^2 - L2.mean[cluster[j],])^2)

  ##Value of target function of Obliclus
  out <- (cr.l + cr.r) / 4

  return(out)
}
