
library("slam")

## test
x <- matrix(c(1, 0, 0, 2, 1, 0), nrow = 3, 
    dimnames = list(A = 1:3, B = 1:2))
s <- as.simple_triplet_matrix(x)
dimnames(x)[[1L]] <- letters[1:3]
names(dimnames(x))[1L] <- 1
x

##
z <- tcrossprod_simple_triplet_matrix(s, x[1:2,])
z

zz <- slam:::.ttcrossprod_simple_triplet_matrix(s, x[1:2,])
identical(z, t(zz))

## bailout
is.na(x) <- 4L

z <- tcrossprod_simple_triplet_matrix(s, x[1:2,])
z

zz <- slam:::.ttcrossprod_simple_triplet_matrix(s, x[1:2,])
identical(z, t(zz))

###
