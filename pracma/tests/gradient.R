##
##  g r a d i e n t . R  Test suite
##


gradient <- pracma::gradient

x <- c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0)
y <- c(1, 2, 3)

Z <- matrix(c(
    1, 1.04, 1.16, 1.36, 1.64,    2,
    4, 4.04, 4.16, 4.36, 4.64,    5,
    9, 9.04, 9.16, 9.36, 9.64,   10), nrow = 3, byrow = TRUE)

X1 <- matrix(c(
    0.04, 0.08, 0.16, 0.24, 0.32, 0.36,
    0.04, 0.08, 0.16, 0.24, 0.32, 0.36,
    0.04, 0.08, 0.16, 0.24, 0.32, 0.36), nrow = 3, byrow = TRUE)

X2 <- matrix(c(
    0.2, 0.4, 0.8, 1.2, 1.6, 1.8,
    0.2, 0.4, 0.8, 1.2, 1.6, 1.8,
    0.2, 0.4, 0.8, 1.2, 1.6, 1.8), nrow = 3, byrow = TRUE)

Y <- matrix(c(
    3, 3, 3, 3, 3, 3,
    4, 4, 4, 4, 4, 4,
    5, 5, 5, 5, 5, 5), nrow = 3, byrow = TRUE)

all.equal(gradient(Z)$X, X1)
all.equal(gradient(Z, x, y)$X, X2)
all.equal(gradient(Z, x, y)$Y, Y)
