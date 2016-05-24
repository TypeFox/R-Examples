
##
library("slam")

s <- as.simple_triplet_matrix(diag(4))
s[1:8] <- 1:8
as.matrix(s)

s[2:3,] <- 1:8
as.matrix(s)

s[,2:3] <- 1:8
as.matrix(s)

s[] <- 1:8
as.matrix(s)

##
local({
    k <- 2:3
    ## Implementing class.
    a <- as.simple_sparse_array(s)
    a[,k]
    a[,k] <- 1:8
    s[,k] <- 1:8
    stopifnot(identical(as.array(a), as.array(s)))
})

###
