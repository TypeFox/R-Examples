
library("geigen")
source("testgv.R")

# tol <-  0 triggers buglet in complex QR
tol <- .Machine$double.eps
set.seed(11)

n <- 5
A.save <- matrix(runif(n*n),nrow=n)
B.save <- matrix(runif(n*n),nrow=n)
A.save <- (A.save+t(A.save))/2
B.save <- (B.save+t(B.save))/2
diag(B.save) <- diag(B.save) + 1

A <- A.save
B <- B.save
B
isSymmetric(A)
isSymmetric(B)
z <- geigen(A,B)
testgv(A,B,z)

A <- A.save+0i
B <- B.save+0i
z <- geigen(A,B)
testgv(A,B,z)

A <- A.save
B <- B.save
z <- geigen(A,B)
testgv(A,B,z)

A <- A.save+0i
B <- B.save+0i
z <- geigen(A,B)
testgv(A,B,z)
