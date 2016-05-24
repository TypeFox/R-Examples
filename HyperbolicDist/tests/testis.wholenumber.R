### Test is.wholenumber
library(HyperbolicDist)

is.wholenumber(-3:5)
is.wholenumber(c(0,0.1,1.3,5))
is.wholenumber(-3:5 + .Machine$double.eps)
is.wholenumber(-3:5 + .Machine$double.eps^0.5)
is.wholenumber(c(2L,3L))
is.wholenumber(c("2L","3L"))
is.wholenumber(0i ^ (-3:3))
is.wholenumber(matrix(1:6, nrow = 3))
is.wholenumber(list(-1:3,2:6))
is.numeric(list(-1:3,2:6))
is.wholenumber(unlist(list(-1:3,2:6)))
