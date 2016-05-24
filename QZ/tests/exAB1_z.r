library(QZ, quiet = TRUE)

### Get eigenvalues
A <- exAB1$A
B <- exAB1$B
ret0 <- qz.zggev(A, B)
ret1 <- qz.zgges(A, B)

# Reordering eigenvalues
S <- ret1$S
T <- ret1$T
Q <- ret1$Q
Z <- ret1$Z
select <- c(TRUE, FALSE, FALSE, TRUE)
ret2 <- qz.ztgsen(S, T, Q, Z, select)

# Verify 0
round(ret0$ALPHA - ret1$ALPHA)
round(ret0$BETA - ret1$BETA)

# Verify 1
A.new <- ret2$Q %*% ret2$S %*% H(ret2$Z)
B.new <- ret2$Q %*% ret2$T %*% H(ret2$Z)
round(A - A.new)
round(B - B.new)

# verify 2
round(ret2$Q %*% H(ret2$Q))
round(ret2$Z %*% H(ret2$Z))
