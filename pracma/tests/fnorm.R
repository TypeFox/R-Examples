##
##  f n o r m  Test suite
##


fnorm <- pracma::fnorm

identical(fnorm(log, sqrt, 1, 2, p = Inf), 1.0)
identical(fnorm(log, sqrt, 1, 2, p = -Inf), sqrt(2) - log(2))
identical(fnorm(log, sqrt, 1, 2, p = 0), Inf)
