##
##  m t r a c e . R  Test suite
##

Trace <- pracma::Trace

identical(Trace(1), 1)
identical(Trace(matrix(c(1,2,3,4,5,6,7,8,9), 3, 3)), 15)
# Error: Trace(matrix(1:12, 3, 4))
