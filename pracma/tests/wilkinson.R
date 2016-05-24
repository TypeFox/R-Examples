###
### wilkinson.R  +++ Test suite +++
###

wilkinson <- pracma::wilkinson

identical(wilkinson(0), NULL)
identical(wilkinson(1), matrix(0, nrow=1, ncol=1))
identical(wilkinson(3), matrix(c(1,1,0, 1,0,1, 0,1,1), 3, 3))
