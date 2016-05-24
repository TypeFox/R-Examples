UK4Apply <-
function(vec, covb, XXSiXi, XSi, Vi, z, n, p)
{
  r1 <- (vec[(n+1):(n+p)] - XSi %*% vec[1:n])
  m <- covb %*% r1
  tlam <- t(vec[1:n] + XXSiXi %*% r1) %*% Vi
  cbind(tlam %*% z,
    sqrt(vec[n+p+1] - tlam %*% vec[1:n] + t(m) %*% vec[(n+1):(n+p)]))
}
