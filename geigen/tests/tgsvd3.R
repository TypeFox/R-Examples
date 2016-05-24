
library("geigen")
source("testgsvd.R")

# double
set.seed(11)
A <- matrix(round(runif(18),2), nrow=3, byrow=TRUE)
B <- matrix(round(runif(30),2), ncol=6, byrow=TRUE)
z <- gsvd(A,B)
testgsvd(z,A,B)

set.seed(11)
A <- matrix(round(runif(18),2), nrow=3, byrow=TRUE)
B <- matrix(round(runif(24),2), ncol=6, byrow=TRUE)
z <- gsvd(A,B)
testgsvd(z,A,B)

set.seed(11)
A <- matrix(round(runif(36),2), ncol=6, byrow=TRUE)
B <- matrix(round(runif(36),2), ncol=6, byrow=TRUE)
z <- gsvd(A,B)
testgsvd(z,A,B)

set.seed(11)
A <- matrix(round(runif(36),2), ncol=6, byrow=TRUE)
B <- matrix(round(runif(24),2), ncol=6, byrow=TRUE)
z <- gsvd(A,B)
testgsvd(z,A,B)

set.seed(11)
A <- matrix(round(runif(12),2), nrow=2, byrow=TRUE)
B <- matrix(round(runif(24),2), ncol=6, byrow=TRUE)
z <- gsvd(A,B)
testgsvd(z,A,B)
