loglikf <-
function(y,w,psi) {
  t(log(t(sapply(y, dnorm, psi$mu,psi$sigma)) %*% psi$p)) %*% w   # %*% is matrix product
}
