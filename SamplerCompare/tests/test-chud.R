# From SamplerCompare, (c) 2010 Madeleine Thompson

# This tests the chud (Cholesky update) and chdd (Cholesky downdate)
# functions in chud.R.

library(SamplerCompare)

### test chud ###

# Define a matrix R and a vector x.  Compute Q so that:
#   t(Q) %*% Q # = t(R) %*% R + x %*% t(x)
# with both chud and by doing a Cholesky factorization.  Fail if
# the two methods do not compute the same answer.

R <- array(c(1, 0, 2, 1), c(2,2))
x <- array(c(0.5, 0.8), c(2,1))


Q <- chud(R,x)
Q.real <- chol(t(R) %*% R + x %*% t(x))

stopifnot(max(abs(Q-Q.real)) < 1e-8)

### test chdd ###

# Use chdd to invert the previous operation on Q.  Fail if we do
# not get the original R back.

R2 <- chdd(Q, x)

stopifnot(max(abs(R-R2)) < 1e-8)
