`orthogonalPolynomials` <-
function(w,x,p) {
  z <- weightedGramSchmidt(outer(x,0:p,"^"),w)$pol[,2:(p+1)]
}

