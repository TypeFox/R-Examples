
library("geigen")
source("testgsvd.R")

# Example from Matlab
# http://nl.mathworks.com/help/matlab/ref/gsvd.html
# also see http://www.netlib.org/lapack/lug/node36.html

A <- matrix(1:15, nrow=5,ncol=3)

B <- matrix(c(8,1,6,
              3,5,7,
              4,9,2), nrow=3)

A
B

z <- gsvd(A,B)
testgsvd(z,A,B)
