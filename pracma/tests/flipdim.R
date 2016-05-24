##
##  f l i p d i m . R  tests
##

flipdim <- pracma::flipdim
flipud <- pracma::flipud
fliplr <- pracma::fliplr
rot90 <- pracma::rot90

a <- matrix(c(1,2,3, 4,5,6, 7,8,9, 10,11,12), nrow=3, ncol=4)
b <- matrix(c(1,2, 3,4), nrow=2, ncol=2, byrow=TRUE)

identical(flipdim(a, 1), flipud(a))
identical(fliplr(a), matrix(c(10,11,12, 7,8,9, 4,5,6, 1,2,3 ), 3, 4))
identical(rot90(b, k=1), matrix(c(2,1, 4,3), 2, 2))
identical(rot90(b, k=6), matrix(c(4,2, 3,1), 2, 2))
identical(rot90(b, k=-1), matrix(c(3,4, 1,2), 2, 2))
