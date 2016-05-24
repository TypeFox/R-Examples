library(QZ, quiet = TRUE)

### http://www.nag.com/lapack-ex/node124.html
A <- exAB1$A
B <- exAB1$B
ret <- qz.zgges(A, B)

# Verify 1
A.new <- ret$Q %*% ret$S %*% H(ret$Z)
B.new <- ret$Q %*% ret$T %*% H(ret$Z)
round(A - A.new)
round(B - B.new)

# verify 2
round(ret$Q %*% H(ret$Q))
round(ret$Z %*% H(ret$Z))
