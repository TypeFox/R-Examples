
library("geigen")
source("testqz.R")

# example from NAG (F08XNF)

A <- matrix(c( -21.1-22.5i , 53.5-50.5i, -34.5+127.5i ,   7.5+.5i,
                 -.46-7.78i, -3.5-37.5i, -15.5+ 58.5i , -10.5-1.5i,
                 4.3-5.5i  , 39.7-17.1i, -68.5+12.5i  ,  -7.5-3.5i,
                 5.5+4.4i  , 14.4+43.3i, -32.5-46.0i  , -19.0-32.5i), nrow=4, byrow=TRUE)

B <- matrix(c( 1.0-5i,  1.6+1.2i , -3.0+0i,       -1i,
               .8-.6i,  3-5i     , -4+3i  , -2.4-3.2i,
               1.0+0i,  2.4+1.8i , -4-5i  ,       -3i,
                   1i, -1.8+2.4i ,   -4i  ,      4-5i ), nrow=4, byrow=TRUE)

z <- gqz(A, B,"R")
testqz(A,B,z)

z <- gqz(A, B)
testqz(A,B,z)

zn <- gqz(A, B, "N")
testqz(A,B,zn)

identical(z,zn)
