
library("geigen")
source("testgsvd.R")

# Example from NAG F08VAF

A <- matrix(c(1,2,3,3,2,1,4,5,6,7,8,8), nrow=4, ncol=3, byrow=TRUE)

B <- matrix(c(-2,-3,3,4,6,5), nrow=2, ncol=3, byrow=TRUE)

A
B

z <- gsvd(A,B)
testgsvd(z,A,B)

# B modified

B <- matrix(1:18,byrow=TRUE, ncol=3)
B

z <- gsvd(A,B)
testgsvd(z,A,B)
