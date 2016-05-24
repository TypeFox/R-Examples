##
library("slam")
x <- matrix(1:12, 4, dimnames = list(NULL, B = 1:3))
s <- as.simple_sparse_array(x)
s

extend_simple_sparse_array(s,  0L)
extend_simple_sparse_array(s, -1L)  ## the same
extend_simple_sparse_array(s,  1L)
extend_simple_sparse_array(s,  2L)
extend_simple_sparse_array(s, -3L)  ## the same

extend_simple_sparse_array(s, c( 0L,  0L))
extend_simple_sparse_array(s, c(-3L, -3L))

## automatic
z <- abind_simple_sparse_array(s, 1:3)
z
all.equal(as.array(z), rbind(x, 1:3))
z <- abind_simple_sparse_array(1:4, s, MARGIN = 2L)
z
all.equal(as.array(z), cbind(1:4, x))

abind_simple_sparse_array(1:3, array(2:4, c(1,3)), array(3:8, c(1,2,3)))
abind_simple_sparse_array(1:3, array(2:4, c(3,1)), array(3:8, c(3,2,1)), MARGIN = 3L) 

## manual
abind_simple_sparse_array(1:3, 2:4)
abind_simple_sparse_array(1:3, 2:4, MARGIN = -1L)
abind_simple_sparse_array(1:3, 2:4, MARGIN = -2L)

###
