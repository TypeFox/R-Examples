analy.deriv <- function(A, T, cluster, info)
{
  N.var <- dim(A)[1]
  N.fac <- dim(A)[2]
  L <- A %*% solve(t(T)) ##Rotated loading matrix
  U <- matrix(0, N.var, info$N.cluster)
  for (i in 1:N.var) {
    U[i, cluster[i]] <- 1
  }

  L2 <- L * L
  one.vec <- rep(1, len=N.fac)

  if (info$method == "oblimin") { ##！！！！！修正が必要！！！！！
    w <- info$oblimin.index
    N <- one.vec %*% t(one.vec) - diag(rep(1, len=N.fac))
    C <- matrix(1/N.var, N.var, N.var)
    Gq <- L * (L2 %*% one.vec %*% t(one.vec) - w * C %*% L2 %*% N -  U %*% solve(t(U) %*% U) %*% t(U) %*% L2)
  }
  else if (info$method == "geomin") {
    delta <- info$geomin.par
    k <- ncol(L)
    p <- nrow(L)
    L2 <- L^2 + delta
    pro <- exp(rowSums(log(L2))/k)
    Gq = (2/k) * (L/L2) * matrix(rep(pro, k), p)  + L * (L2 - U %*% solve(t(U) %*% U) %*% t(U) %*% L2)
  }

  return(- t(t(L) %*% Gq %*% solve(T)))
}
