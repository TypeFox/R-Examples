#  File tests/valued_fit.R in package ergm.count, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2008-2016 Statnet Commons
#######################################################################
library(ergm.count)
set.seed(0)

## Poisson-reference ERGM

n <- 5

m <- matrix(rpois(n^2,2),n,n)
diag(m) <- 0
y <- as.network(m, matrix.type="a", directed=TRUE, ignore.eval=FALSE, names.eval="w")

truth <- log(sum(m)/n/(n-1))
diag(m) <- NA

efit <- ergm(y ~ sum, response="w", reference=~Poisson, verbose=TRUE)

summary(efit)

true.llk <- sum(dpois(na.omit(c(m)), exp(coef(efit)), log=TRUE)) - sum(dpois(na.omit(c(m)), 1, log=TRUE))

stopifnot(abs(coef(efit)-truth)<0.01)
stopifnot(abs(true.llk - logLik(efit))<0.2)
