##
library("slam")
## sparse
x <- simple_sparse_zero_array(dim = c(3, 4, 2))
				    ## removal of subscripts
k <- matrix(c(2, 1, 1, 0, 1, 1), c(2, 3), byrow = TRUE)
k
x[k] <- 1
x[3, 1, 1] <- 2 
x[c(17, 17)] <- c(2, 3)		    ## duplicate subscripts
x[c(5, NA, 6)] <- 4		    ## recycling
data.frame(v = x$v, i = x$i,
    k = .Call("_vector_index", x$dim, x$i))

##
x[, -1, 1] <- 0			    ## zero elements
data.frame(v = x$v, i = x$i,
    k = .Call("_vector_index", x$dim, x$i))
x[-c(2, 3)] <- 0
data.frame(v = x$v, i = x$i,
    k = .Call("_vector_index", x$dim, x$i))

x[] <- 0
str(x)

## misc
x[integer()] <- 1
x[matrix(integer(), nrow = 0, ncol = 3)] <- 1
str(x)

try(x[c(NA, 2, 3)] <- 1:2)	    ## not allowed
## works with R >= 3.x
try(x[-c(.Machine$integer.max + 1, 1)] <- c(1, 2))


## reference
x <- matrix(1:6, nrow = 3)

## matrix indexing
k <- matrix(c(1, 1, 2, 2, 1, 1), ncol = 2, byrow = TRUE)
k

z <- x
z[k] <- -1
z

z <- x
z[k] <- -(1:3)			    ## last in sequence 
z

## implicit vector indexing
k <- matrix(k, nrow = 2)
as.vector(k)

z <- x
z[k] <- -1
z

z <- x
z[k] <- -(1:6)			    ## last in sequence
z

## missing values
z <- x
z[c(NA, 1, 2)] <- -1
z

z <- x
try(z[c(NA, 1, 2)] <- -(1:2))	    ## not allowed

k[1L] <- NA			    ## implicit vector indexing
as.vector(k)
z <- x
z[k] <- -1
z

k <- matrix(c(NA, 1, 1, 1, 2, 2), ncol = 2, byrow = TRUE)
k

z <- x
z[k] <- -1
z

z <- x
try(z[k] <- -(1:2))		    ## not allowed

## zeros
z <- x
z[c(0, 1)] <- -1
z

z <- x
z[c(0, 1)] <- -(1:2)
z

k <- matrix(c(1, 1, 0, 2), ncol = 2, byrow = TRUE)
k

z <- x
z[k] <- -1
z

z <- x
z[k] <- -(1:2)
z

## extending
k <- matrix(c(1, 4), ncol = 2)

z <- x
try(z[k] <- 1)			    ## not allowed

z[c(1, 8)] <- 1			    ## not implemented
z

## misc
z <- x
try(z[-c(.Machine$integer.max + 1, 1)] <- c(1, 2))

###
