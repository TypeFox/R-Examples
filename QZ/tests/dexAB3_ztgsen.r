library(QZ, quiet = TRUE)

### http://www.nag.com/numeric/fl/nagdoc_fl23/xhtml/F08/f08yuf.xml
S <- exAB3$S
T <- exAB3$T
Q <- exAB3$Q
Z <- exAB3$Z
select <- c(FALSE, TRUE, TRUE, FALSE)
ret <- qz.ztgsen(S, T, Q, Z, select)

# Verify 1
S.new <- ret$Q %*% ret$S %*% H(ret$Z)
T.new <- ret$Q %*% ret$T %*% H(ret$Z)
round(S - S.new)
round(T - T.new)

# verify 2
round(ret$Q %*% H(ret$Q))
round(ret$Z %*% H(ret$Z))
