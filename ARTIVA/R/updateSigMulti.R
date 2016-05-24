updateSigMulti <-
function(phase, X, Y, E, Sall, Ball, Sig2, Mphase, alphad2, betad2, v0, gamma0){
  
  posPhase = which(E == phase)
  S = Sall[posPhase,]
  k = sum(Sall[posPhase,])-1
  y = Y[(Mphase[phase]:(Mphase[E[posPhase+1]]-1))]
  x = X[(Mphase[phase]:(Mphase[E[posPhase+1]]-1)),]
  delta2 = rinvgamma(1, shape=k + alphad2, scale=betad2 + Ball[posPhase, which(S == 1)] %*% t(x[, which(S == 1)]) %*% x[,which(S == 1)] %*% Ball[posPhase,which(S == 1)] / (2 * Sig2) )
  matPx = computePx(length(y), x[, which(S == 1)], delta2)
    
  total = t(y) %*% matPx %*%y
  out = rinvgamma(1, shape=v0/2 + length(y)/2, scale=(gamma0 + total)/2)
  return(out )
  
}
