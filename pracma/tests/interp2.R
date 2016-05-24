##
##  i n t e r p 2 . R Test suite
##


interp2 <- pracma::interp2

x <- c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0)
y <- c(1, 2, 3)

Z <- matrix(c(
    1, 1.04, 1.16, 1.36, 1.64,    2,
    4, 4.04, 4.16, 4.36, 4.64,    5,
    9, 9.04, 9.16, 9.36, 9.64,   10), nrow = 3, byrow = TRUE)

all.equal(interp2(x, y, Z, 0.55, 2.55, method = "constant"), 4.16)
all.equal(interp2(x, y, Z, 0.55, 2.55, method = "nearest"),  9.36)

all.equal(interp2(x, y, Z, 0.5,  2.5,  method = "linear"), 6.76)
all.equal(interp2(x, y, Z, 0.55, 2.5,  method = "linear"), 6.81)
all.equal(interp2(x, y, Z, 0.5,  2.55, method = "linear"), 7.01)
all.equal(interp2(x, y, Z, 0.55, 2.55, method = "linear"), 7.06)

all.equal(interp2(x, y, Z, 0.0, 1.5, method = "linear"), 2.5)
all.equal(interp2(x, y, Z, 0.1, 1.0, method = "linear"), 1.02)
all.equal(interp2(x, y, Z, 1.0, 2.5, method = "linear"), 7.5)
all.equal(interp2(x, y, Z, 0.9, 3.0, method = "linear"), 9.82)
