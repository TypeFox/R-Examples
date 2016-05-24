library(QZ, quiet = TRUE)

### Get eigenvalues
A <- exAB2$A
B <- exAB2$B
ret0 <- qz.dggev(A, B)
ret1 <- qz.dgges(A, B)

# Reordering eigenvalues
S <- ret1$S
T <- ret1$T
Q <- ret1$Q
Z <- ret1$Z
select <- c(FALSE, TRUE, TRUE, FALSE)
ret2 <- qz.dtgsen(S, T, Q, Z, select)

# Verify 0
round(ret0$ALPHAR - ret1$ALPHAR)
round(ret0$ALPHAI - ret1$ALPHAI)
round(ret0$BETA - ret1$BETA)

# Verify 1
A.new <- ret2$Q %*% ret2$S %*% t(ret2$Z)
B.new <- ret2$Q %*% ret2$T %*% t(ret2$Z)
round(A - A.new)
round(B - B.new)

# verify 2
round(ret2$Q %*% t(ret2$Q))
round(ret2$Z %*% t(ret2$Z))
