##
## https://stat.ethz.ch/pipermail/r-help/2007-September/141324.html
##

library(quadprog)

for(n in 1:10){

  Dmat <- diag(2*n)
  dvec <- rep(0, 2*n)
  
  Amat <- cbind(rep(c(1,-1), n),
                rep(1, 2*n),
                kronecker(diag(n), diag(c(-1,1))))
  bvec <- c(1, 1, rep(c(-1,0), n))
  
  meq <- 2

  cat("\n\nFor problem of size", n,":\n")
  print(lapply(solve.QP(Dmat, dvec, Amat, bvec, meq), zapsmall))
}
