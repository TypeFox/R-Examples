sampleDelta2 <-
function(pos, x, q, B, S, sig2, alphad2, betad2){
  # INPUT: pos, the considered state
  #        xPos, the observations of X in state i
  #        B,S,Sig2
  #        alphad2,betad2.
  # OUTPUT: delta2
  # depends on: . 
  plus = 0
  if(sum(S[pos,])>0){
    Bi = B[pos, which(S[pos,] == 1)]
    xi = x[, which(S[pos,] == 1)]
    plus = Bi %*% t(xi) %*% xi %*% Bi / (2* sig2)
  }
  out = rinvgamma(1, shape=sum(S[pos,1:q]) + alphad2, scale=betad2 + plus)
  
  return(out)
}
