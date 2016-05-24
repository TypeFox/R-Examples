
library("geigen")
source("testgv.R")

# tol <-  0 triggers buglet in complex QR
tol <- .Machine$double.eps
set.seed(11)

n <- 5
A.save <- matrix(runif(n*n),nrow=n)
B.save <- matrix(runif(n*n),nrow=n)

A <- A.save
B <- B.save
B[5,] <- (B[4,]+B[3,])/2 + tol
z <- geigen(A,B)
testgv(A,B,z)

A <- A.save+0i
B <- B.save+0i
B[5,] <- (B[4,]+B[3,])/2 + tol
z <- geigen(A,B)
testgv(A,B,z)

A <- A.save
B <- B.save
B[1,] <- (B[2,]+B[3,])/2 + tol
z <- geigen(A,B)
testgv(A,B,z)

A <- A.save+0i
B <- B.save+0i
B[1,] <- (B[2,]+B[3,])/2 + tol
z <- geigen(A,B)
testgv(A,B,z)
