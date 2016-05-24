library(QZ, quiet = TRUE)

### http://www.nag.com/lapack-ex/node92.html
A <- exA1$A
ret <- qz.zgeev(A)

# Verify 1
diff.R <- A %*% ret$V - matrix(ret$W, 4, 4, byrow = TRUE) * ret$V
diff.L <- H(ret$U) %*% A - matrix(ret$W, 4, 4) * H(ret$U)
round(diff.R)
round(diff.L)

# Verify 2
round(ret$U %*% H(ret$U))
round(ret$V %*% H(ret$V))
