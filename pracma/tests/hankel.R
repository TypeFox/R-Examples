##
##  hankel.R  Test
##

hankel <- pracma::hankel

identical(hankel(2), matrix(2, nrow=1, ncol=1))
identical(hankel(1:3), matrix(c(1,2,3,2,3,0,3,0,0), 3, 3))
identical(hankel(1:3, 3:1), matrix(c(1,2,3,2,3,2,3,2,1), 3, 3))
identical(hankel(1:3, 2:1), matrix(c(1,2,3,2,3,1), 3, 2))
identical(hankel(1:2, 3:1), matrix(c(1,2,2,2,2,1), 2, 3))
