
## test C interfaces

library(proxy)

set.seed(20070630)

x <- as.dist(matrix(runif(25),5,5))
x
attributes(x)

z <- .Call("R_subset_dist", x, 3)
z

unclass(z)

.Call("R_subset_dist", x, c(1,3,5))

attr(x, "Labels") <- LETTERS[1:5]

z <- .Call("R_subset_dist", x, c("A","C","E"))
z
attributes(z)

attr(x, "Labels") <- NULL

.Call("R_rowSums_dist", x, FALSE)
.Call("R_rowSums_dist", z, FALSE)

.Call("R_row_dist", x, FALSE)       # row()
.Call("R_row_dist", x, TRUE)        # col()

## test R interfaces

dim(x)
dimnames(x) <- letters[1:5]
dimnames(x)
names(x)    <- LETTERS[1:5]
names(x)

row.dist(x)
col.dist(x)

subset(x, c(1,3,5))
x[[c(1,3,5)]]
x[c(1,3,5)]                         # as usual

x[[-1]]                             # drop subscripts

x[[1]]                              # empty

###
