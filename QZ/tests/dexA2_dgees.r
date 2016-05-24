library(QZ, quiet = TRUE)

### http://www.nag.com/lapack-ex/node89.html
A <- exA2$A
ret <- qz.dgees(A)

# Verify 1
A.new <- ret$Q %*% ret$T %*% solve(ret$Q)
round(A - A.new)

# verify 2
round(ret$Q %*% solve(ret$Q))
