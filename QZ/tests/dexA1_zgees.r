library(QZ, quiet = TRUE)

### http://www.nag.com/lapack-ex/node94.html
A <- exA1$A
ret <- qz.zgees(A)

# Verify 1
A.new <- ret$Q %*% ret$T %*% H(ret$Q)
round(A - A.new)

# verify 2
round(ret$Q %*% H(ret$Q))
