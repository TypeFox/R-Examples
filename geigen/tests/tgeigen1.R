
# from Octave's eigen test
library(geigen)
source("testgv.R")

A <- matrix(c(1, 2, -1, 1),nrow=2, byrow=TRUE)
B <- matrix(c(3, 3,  1, 2),nrow=2, byrow=TRUE)

z <- geigen(A,B)
testgv(A,B,z)

A <- matrix(c(1, 2,  2, 1),nrow=2, byrow=TRUE)
B <- matrix(c(3,-2, -2, 3),nrow=2, byrow=TRUE)

z <- geigen(A,B)
testgv(A,B,z)

z <- geigen(A,B, symmetric=TRUE)
testgv(A,B,z)

# Complex
A <- matrix(c(1+3i, 2+1i, 2-1i, 1+3i),nrow=2, byrow=TRUE)
B <- matrix(c(5+9i, 2+1i, 2-1i, 5+9i),nrow=2, byrow=TRUE)

z <- geigen(A,B)
testgv(A,B,z)

# Hermitian 
A <- matrix(c(3, 2+1i, 2-1i, 5),nrow=2, byrow=TRUE)
B <- matrix(c(5, 2+1i, 2-1i, 5),nrow=2, byrow=TRUE)

z <- geigen(A,B)
testgv(A,B,z)

z <- geigen(A,B, symmetric=TRUE)
testgv(A,B,z)
 
A <- matrix(c(1+3i, 2+3i, 3-8i, 8+3i),nrow=2, byrow=TRUE)
B <- matrix(c(8+1i, 3+1i, 4-9i, 3+1i ),nrow=2, byrow=TRUE)

z <- geigen(A,B)
testgv(A,B,z)

A <- matrix(c(1, 2, 3, 8),nrow=2, byrow=TRUE)
B <- matrix(c(8, 3, 4, 3),nrow=2, byrow=TRUE)

z <- geigen(A,B)
testgv(A,B,z)
