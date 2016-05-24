library(QZ, quiet = TRUE)

### Get eigenvalues
A <- exA1$A
ret0 <- qz.zgeev(A)
ret1 <- qz.zgees(A)

# Reordering eigenvalues
T <- ret1$T
Q <- ret1$Q
select <- c(TRUE, FALSE, FALSE, TRUE)
ret2 <- qz.ztrsen(T, Q, select)

# Verify 0
round(ret0$W - ret1$W)

# Verify 1
A.new <- ret2$Q %*% ret2$T %*% solve(ret2$Q)
round(A - A.new)

# verify 2
round(ret2$Q %*% solve(ret2$Q))
