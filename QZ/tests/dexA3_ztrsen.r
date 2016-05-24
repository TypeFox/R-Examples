library(QZ, quiet = TRUE)

### http://www.nag.com/numeric/fl/nagdoc_fl23/xhtml/F08/f08quf.xml
T <- exA3$T
Q <- exA3$Q
select <- c(TRUE, FALSE, FALSE, TRUE)
ret <- qz.ztrsen(T, Q, select)

# Verify 1
A <- Q %*% T %*% solve(Q)
A.new <- ret$Q %*% ret$T %*% solve(ret$Q)
round(A - A.new)

# verify 2
round(ret$Q %*% solve(ret$Q))
