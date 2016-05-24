##
##  c h e b P o l y . R  Test suite
##


chebPoly <- pracma::chebPoly
chebCoeff <- pracma::chebCoeff
chebApprox <- pracma::chebApprox

identical(chebPoly(6),
          matrix( c(0,   0,   0,   0,   0,   0,   1,
                    0,   0,   0,   0,   0,   1,   0,
                    0,   0,   0,   0,   2,   0,  -1,
                    0,   0,   0,   4,   0,  -3,   0,
                    0,   0,   8,   0,  -8,   0,   1,
                    0,  16,   0, -20,   0,   5,   0,
                   32,   0, -48,   0,  18,   0,  -1),
                 nrow = 7, ncol = 7, byrow=TRUE))

f <- function(x) 1 + x/1 + x^2/2 + x^3/6 + x^4/24 + x^5/120 + x^6/720
cC <- chebCoeff(f, -1, 1, 6)
cC[1] <- cC[1]/2
all.equal(cC,
          c(1.26606, 1.13021, 0.27148, 0.04427, 0.00547, 0.00052, 0.00004),
          tol = 1e-5)

x <- seq(-1, 1, length.out=7)
y <- chebApprox(x, function(x) x^2, -1, 1, 6)
all.equal(x^2, y, tol = 1e-7)
