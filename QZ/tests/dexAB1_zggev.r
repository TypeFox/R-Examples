library(QZ, quiet = TRUE)

### http://www.nag.com/lapack-ex/node122.html
A <- exAB1$A
B <- exAB1$B
ret <- qz.zggev(A, B)

# Verify 1
(lambda <- ret$ALPHA / ret$BETA)    # Unstable
diff.R <- matrix(ret$BETA, 4, 4, byrow = TRUE) * A %*% ret$V -
          matrix(ret$ALPHA, 4, 4, byrow = TRUE) * B %*% ret$V
diff.L <- matrix(ret$BETA, 4, 4) * H(ret$U) %*% A -
          matrix(ret$ALPHA, 4, 4) * H(ret$U) %*% B
round(diff.R)
round(diff.L)

# Verify 2
round(ret$U %*% solve(ret$U))
round(ret$V %*% solve(ret$V))
