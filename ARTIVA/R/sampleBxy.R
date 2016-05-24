sampleBxy <-
function(xi, y, Sig2, delta2){
  # INPUT: xi, yi, sig2i, delta2
  # OUTPUT: B
  # depends on: .
  Ml = (delta2 / (delta2+1)) * pseudoinverse(t(xi) %*% xi)
  out = mvrnorm(1, mu=Ml %*% t(xi) %*% y, Sigma=Sig2*Ml)
  return(out)
}
