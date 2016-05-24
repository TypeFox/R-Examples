getNormalizer <-
function(avgmat, perm){
  factorA = sum(diag((avgmat$B-avgmat$G) %*% avgmat$B)) /(avgmat$m-1)
  P = perm2Mat(perm)
  factorB = sum(diag((avgmat$B-avgmat$G) %*% P %*% avgmat$B %*% t(P))) /(avgmat$m-1)
  return(list(facA = factorA, facB = factorB, norm = 1/(factorA-factorB)))
}
