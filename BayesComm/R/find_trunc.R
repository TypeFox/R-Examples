find_trunc <-
function (mu, y) {
  n <- dim(mu)[1]
  nsp <- dim(mu)[2]
  eps <- 10 ^ -7
  trunc <- array(-Inf, dim = c(n, nsp, 2))
  trunc[, , 2] <- -mu - eps
  ind <- which(t(y) == 1)
  trunc[, , 1][ind] <- (eps - mu)[ind]
  trunc[, , 2][ind] <- Inf
  trunc
}
