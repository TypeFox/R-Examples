
library(potts)

set.seed(42)

ncolor <- as.integer(2)
nrow <- sample(100:200, 1)
ncol <- sample(100:200, 1)
x <- matrix(sample(ncolor, nrow * ncol, replace = TRUE), nrow = nrow)
foo <- packPotts(x, ncolor)
foo[2] == 3 # log2pixelsperbyte
bar <- inspectPotts(foo)
identical(nrow, bar$nrow)
identical(ncol, bar$ncol)
identical(ncolor, bar$ncolor)
baz <- unpackPotts(foo)
identical(x, baz)

ncolor <- as.integer(3)
nrow <- sample(100:200, 1)
ncol <- sample(100:200, 1)
x <- matrix(sample(ncolor, nrow * ncol, replace = TRUE), nrow = nrow)
foo <- packPotts(x, ncolor)
foo[2] == 2 # log2pixelsperbyte
bar <- inspectPotts(foo)
identical(nrow, bar$nrow)
identical(ncol, bar$ncol)
identical(ncolor, bar$ncolor)
baz <- unpackPotts(foo)
identical(x, baz)

ncolor <- as.integer(4)
nrow <- sample(100:200, 1)
ncol <- sample(100:200, 1)
x <- matrix(sample(ncolor, nrow * ncol, replace = TRUE), nrow = nrow)
foo <- packPotts(x, ncolor)
foo[2] == 2 # log2pixelsperbyte
bar <- inspectPotts(foo)
identical(nrow, bar$nrow)
identical(ncol, bar$ncol)
identical(ncolor, bar$ncolor)
baz <- unpackPotts(foo)
identical(x, baz)

ncolor <- as.integer(5)
nrow <- sample(100:200, 1)
ncol <- sample(100:200, 1)
x <- matrix(sample(ncolor, nrow * ncol, replace = TRUE), nrow = nrow)
foo <- packPotts(x, ncolor)
foo[2] == 1 # log2pixelsperbyte
bar <- inspectPotts(foo)
identical(nrow, bar$nrow)
identical(ncol, bar$ncol)
identical(ncolor, bar$ncolor)
baz <- unpackPotts(foo)
identical(x, baz)

ncolor <- as.integer(17)
nrow <- sample(100:200, 1)
ncol <- sample(100:200, 1)
x <- matrix(sample(ncolor, nrow * ncol, replace = TRUE), nrow = nrow)
foo <- packPotts(x, ncolor)
foo[2] == 0 # log2pixelsperbyte
bar <- inspectPotts(foo)
identical(nrow, bar$nrow)
identical(ncol, bar$ncol)
identical(ncolor, bar$ncolor)
baz <- unpackPotts(foo)
identical(x, baz)

ncolor <- as.integer(256)
nrow <- sample(100:200, 1)
ncol <- sample(100:200, 1)
x <- matrix(sample(ncolor, nrow * ncol, replace = TRUE), nrow = nrow)
foo <- packPotts(x, ncolor)
foo[2] == 0 # log2pixelsperbyte
bar <- inspectPotts(foo)
identical(nrow, bar$nrow)
identical(ncol, bar$ncol)
identical(ncolor, bar$ncolor)
baz <- unpackPotts(foo)
identical(x, baz)

ncolor <- as.integer(257)
nrow <- sample(100:200, 1)
ncol <- sample(100:200, 1)
x <- matrix(sample(ncolor, nrow * ncol, replace = TRUE), nrow = nrow)
try(foo <- packPotts(x, ncolor)) # should fail

