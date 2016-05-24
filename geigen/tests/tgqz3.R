
library("geigen")
source("testqz.R")

# example from NAG

A <- matrix(c(  3.9, 12.5,-34.5,-0.5,
                4.3, 21.5,-47.5, 7.5,
                4.3, 21.5,-43.5, 3.5,
                4.4, 26.0,-46.0, 6.0), nrow=4, byrow=TRUE)

B <- matrix(c( 1.0, 2.0, -3.0, 1.0,
               1.0, 3.0, -5.0, 4.0,
               1.0, 3.0, -4.0, 3.0,
               1.0, 3.0, -4.0, 4.0), nrow=4, byrow=TRUE)

z <- gqz(A, B,"R")
testqz(A,B,z)
