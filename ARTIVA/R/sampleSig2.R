sampleSig2 <-
function(y, Px, v0, gamma0) {
  out = rinvgamma(1, shape=v0/2 + length(y)/2, scale = (gamma0 + t(y) %*% Px %*% y)/2)
  return(out)
}
