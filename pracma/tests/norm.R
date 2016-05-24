##
##  n o r m  Test suite
##


Norm <- pracma::Norm

identical(Norm(c(3, 4)), 5)
identical(Norm(c(1, 1, 1), p=2), sqrt(3))
identical(Norm(1:10, p = 1), sum(1:10)+0.0)
identical(Norm(1:10, p = 0), Inf)
identical(Norm(1:10, p = Inf), max(1:10))
identical(Norm(1:10, p = -Inf), min(1:10))
