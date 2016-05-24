##
##  q u a d . R  Test suite
##


quad <- pracma::quad

all.equal(quad(sin, 0, pi), 2, tol = 1e-7)
all.equal(quad(sin, 0, 2*pi), 0, tol = 1e-7)
all.equal(quad(exp, 0, 1), exp(1) - 1, tol=1e-7)
