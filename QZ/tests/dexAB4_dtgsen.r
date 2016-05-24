library(QZ, quiet = TRUE)

### http://www.nag.com/numeric/fl/nagdoc_fl23/xhtml/F08/f08ygf.xml
S <- exAB4$S
T <- exAB4$T
Q <- exAB4$Q
Z <- exAB4$Z
select <- c(TRUE, FALSE, FALSE, TRUE)
ret <- qz.dtgsen(S, T, Q, Z, select)

# Verify 1
S.new <- ret$Q %*% ret$S %*% t(ret$Z)
T.new <- ret$Q %*% ret$T %*% t(ret$Z)
round(S - S.new)
round(T - T.new)

# verify 2
round(ret$Q %*% t(ret$Q))
round(ret$Z %*% t(ret$Z))
