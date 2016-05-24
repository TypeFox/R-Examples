## Test results before and after loading gdata

m <- factor(c('a','b','c'))

( m1 <- reorder(m, X=c(3, 2, 1)) )

library(gdata)

( m2 <- reorder(m, X=c(3, 2, 1)) )

stopifnot(identical(m1,m2))
