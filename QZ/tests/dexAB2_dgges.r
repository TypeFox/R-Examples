library(QZ, quiet = TRUE)

### http://www.nag.com/lapack-ex/node119.html
A <- exAB2$A
B <- exAB2$B
ret <- qz.dgges(A, B)

# Verify 1
A.new <- ret$Q %*% ret$S %*% t(ret$Z)
B.new <- ret$Q %*% ret$T %*% t(ret$Z)
round(A - A.new)
round(B - B.new)

# verify 2
round(ret$Q %*% t(ret$Q))
round(ret$Z %*% t(ret$Z))
