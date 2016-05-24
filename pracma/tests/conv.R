##
##  c o  n v . r  Test suite
##


conv <- pracma::conv
deconv <- pracma::deconv

all.equal(conv(c(1, 1, 1), 1), c(1, 1, 1))
all.equal(conv(c(1, 1, 1), c(0, 0, 1)), c(0, 0, 1, 1, 1))
all.equal(conv(c(-0.5, 1, -1), c(0.5, 0, 1)), c(-0.25, 0.5, -1, 1, -1))

b <- c(-0.25, 0.5, -1, 1, -1)
a <- c(0.5, 0, 1)
d <- deconv(b, a)
all.equal(d$q, c(-0.5, 1, -1))
all.equal(d$r, c(0, 0))
