
##
library(slam)

## zero dimension
a <- as.simple_sparse_array(array(0L, 0L))
drop_simple_sparse_array(a)

## invalid
a <- simple_sparse_array(rep(1L, 2L), c(1L, -1L))
a <- reduce_simple_sparse_array(a)
as.array(a)

## not minimal
x <- matrix(1:6, 3L, 2, dimnames = list(NULL, NULL))
a <- as.simple_sparse_array(x)
z <- reduce_simple_sparse_array(a)
identical(a, z)

##
v <- c("logical", "integer", "double", "complex", "character", "list")
stopifnot(any(sapply(v, function(v) 
    !.Call("__valid_v", vector(typeof(v), 1L)))))

##
