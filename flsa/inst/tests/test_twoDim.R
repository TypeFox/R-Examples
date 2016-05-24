library(flsa)
library(testthat)

res.twodim <- flsa(matrix(1:16, nrow=4))
res.twodim.sol <- flsaGetSolution(res.twodim, lambda2=c(0, 1, 2, 4, 8))

## check that the calculated lambda values are correct
expect_equal(unique(res.twodim$BeginLambda), c(0, 1, 2, 4, 8))
expect_equal(as.vector(res.twodim.sol[1, ]), as.double(1:16), tolerance=0.0001)
expect_equal(as.vector(res.twodim.sol[2, ]), rep(c(3.0, 4.0, 6.0, 7.0, 10.0, 11.0, 13.0, 14.0), each=2), tolerance=0.00001)
expect_equal(as.vector(res.twodim.sol[3, ]), rep(c(4.5, 6.5, 10.5, 12.5), each=4), tolerance=0.000001)
expect_equal(as.vector(res.twodim.sol[4, ]), rep(c(6.5, 10.5), each=8), tolerance=0.000001)
expect_equal(as.vector(res.twodim.sol[5, ]), rep(8.5, each=16), tolerance=0.0000001)
