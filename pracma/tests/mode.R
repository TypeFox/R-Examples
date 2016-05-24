##
##  m o d e . R
##


Mode <- pracma::Mode

x <- c(1:100, rep(5,3), rep(27,5), rep(71,4), rep(89,2), rep(100, 5))
identical(Mode(x), 27)

x <- as.factor(x)
identical(Mode(x), "27")
