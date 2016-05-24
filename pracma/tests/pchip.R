##
##  p c h i p . R  Test suite
##


pchip <- pracma::pchip

x <- c(1, 2, 3, 4, 5, 6)
y <- c(16, 18, 21, 17, 15, 12)

xs <- c(1.5, 2.5, 3.5, 4.5, 5.5)
ys <- pchip(x, y, xs)
# ys <- interp1(x, y, xs, method="cubic")       # the same
# 16.88750 19.80000 19.33333 15.96667 13.63750

yml <- c(16.887499999999999, 19.800000000000001, 19.333333333333332,
         15.966666666666667, 13.637499999999999)

all.equal(ys, yml, tolerance = 1e-7)
