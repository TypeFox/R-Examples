
library(geigen)
source("testgv.R")

# Data from R-help mailinglist: Solve an ordinary or generalized eigenvalue problem in R?
# started on 19-04-2012

A <- matrix(c(1457.738, 1053.181, 1256.953,
              1053.181, 1213.728, 1302.838,
              1256.953, 1302.838, 1428.269), nrow=3, byrow=TRUE)

B <- matrix(c(4806.033, 1767.480, 2622.744,
              1767.480, 3353.603, 3259.680,
              2622.744, 3259.680, 3476.790), nrow=3, byrow=TRUE)

z1 <- geigen(A, B, symmetric=FALSE, only.values=TRUE)
z2 <- geigen(A, B, symmetric=FALSE, only.values=FALSE )
all.equal(z1$values, z2$values)
testgv(A,B,z2)  

# geigen(A, B)
z1 <- geigen(A, B, only.values=TRUE)
z2 <- geigen(A, B, only.values=FALSE)
all.equal(z1$values, z2$values)
testgv(A,B,z2)  

A.c <- A + 1i
B.c <- B + 1i
isSymmetric(A.c)
isSymmetric(B.c)

z1 <- geigen(A.c, B.c, only.values=TRUE)
z2 <- geigen(A.c, B.c, only.values=FALSE)
all.equal(z1$values, z2$values)
testgv(A.c,B.c,z2)

A[upper.tri(A)] <- A[upper.tri(A)] + 1i   
A[lower.tri(A)] <- Conj(t(A[upper.tri(A)]))

B[upper.tri(B)] <- B[upper.tri(B)] + 1i
B[lower.tri(B)] <- Conj(t(B[upper.tri(B)])) 

isSymmetric(A)
isSymmetric(B) 

z1 <- geigen(A, B, only.values=TRUE)
z2 <- geigen(A, B, only.values=FALSE)
all.equal(z1$values, z2$values)
testgv(A,B,z2)  
