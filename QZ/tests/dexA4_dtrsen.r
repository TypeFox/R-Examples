library(QZ, quiet = TRUE)

### http://www.nag.com/numeric/fl/nagdoc_fl22/xhtml/F08/f08qgf.xml
T <- exA4$T
Q <- exA4$Q
select <- c(TRUE, FALSE, FALSE, TRUE)
ret <- qz.dtrsen(T, Q, select)

# Verify 1
A <- Q %*% T %*% solve(Q)
A.new <- ret$Q %*% ret$T %*% solve(ret$Q)
round(A - A.new)

# verify 2
round(ret$Q %*% t(ret$Q))
