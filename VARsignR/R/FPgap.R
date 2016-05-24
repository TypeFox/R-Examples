FPgap <-
function(g, rescale, targets, impulses, nstep){
  gap <- 0.0
  a <- matrix(stereo(v=g))
  for(k in 1:nstep){
    scal  <-  matrix(rescale[k,] * (targets[k,] - impulses[k, , ] %*% a))
    scal <- sum(scal^2)
    gap <- gap + scal
  }
  return(gap)
}
