pointfit <- function(xi,x) {  # row vectors
  xiq <- apply(xi,2,mean)
  xq  <- apply(x ,2,mean)
  A   <- sweep(x ,2,xq )
  B   <- sweep(xi,2,xiq)
  sv  <- La.svd(t(A) %*% B)
  Q   <- t(sv$vt) %*% t(sv$u)
#  f   <- mean((B/(A %*% t(Q))))
  lf  <- log(B/(A %*% t(Q)))
  lf  <- lf[is.finite(lf)]
  f   <- exp(mean(lf))
  Qf   <- Q*f
  list(Q = Qf, tr = as.vector(xiq - xq %*% Qf), factor = f, res = B - A %*% t(Qf))
}
