###
###     h y p o t  Tests
###


hypot <- pracma::hypot
identical(hypot(3,4), 5)
identical(hypot(c(0,0), c(3,4)), c(3,4))
