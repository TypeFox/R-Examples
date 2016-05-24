library(QZ, quiet = TRUE)

### http://www.nag.com/lapack-ex/node87.html
A <- exA2$A
ret <- qz.dgeev(A)

# Verify 1
diff.R <- A %*% ret$V - matrix(ret$W, 4, 4, byrow = TRUE) * ret$V
diff.L <- t(ret$U) %*% A - matrix(ret$W, 4, 4) * t(ret$U)
round(diff.R)
round(diff.L)

# Verify 2
round(ret$U %*% solve(ret$U))
round(ret$V %*% solve(ret$V))
