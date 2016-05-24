##
##  c o  n d . r  Test suite
##


cond <- pracma::cond
normest <- pracma::normest

hilb <- pracma::hilb
all.equal(c(cond(hilb(1)), cond(hilb(2)), cond(hilb(3)), cond(hilb(4))),
          c(1, 19.281470, 524.056778,  15513.738739),
          tolerance = 1e-6)
magic <- pracma::magic
all.equal(normest(magic(5)), max(svd(magic(5))$d))
all.equal(normest(pracma::magic(100)), 500050)
