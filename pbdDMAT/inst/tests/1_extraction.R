library(pbdTEST)
settings(mpi=TRUE)

.BLDIM <- 2
seed <- 10


### --------------------------------------
module("Extraction: 1 column")

x <- matrix(1:20, ncol=1)
dx <- as.ddmatrix(x)

test("[1:6, ]", {
  a <- x[1:6, , drop=FALSE]
  b <- as.matrix(dx[1:6, ])
})

collect()



### --------------------------------------
module("Extraction: 1 row")

x <- matrix(1:20, nrow=1)
dx <- as.ddmatrix(x)

test("[, 1:6]", {
  a <- x[, 1:6, drop=FALSE]
  b <- as.matrix(dx[, 1:6])
})

collect()



### --------------------------------------
module("Extraction: square")

x <- matrix(1:100, ncol=10)
dx <- as.ddmatrix(x, 2)

test("[1:3, 1:2]", {
  a <- x[1:3, 1:2, drop=FALSE]
  b <- as.matrix(dx[1:3, 1:2])
})

test("[1, 1]", {
  a <- x[1, 1, drop=FALSE]
  b <- as.matrix(dx[1, 1])
})

test("[7, 7]", {
  a <- x[7, 7, drop=FALSE]
  b <- as.matrix(dx[7, 7])
})

test("[1:2, 1:2]", {
  a <- x[1:2, 1:2, drop=FALSE]
  b <- as.matrix(dx[1:2, 1:2])
})

test("[-1, -1]", {
  a <- x[-1, -1, drop=FALSE]
  b <- as.matrix(dx[-1, -1])
})

test("[5:7, 8:9]", {
  a <- x[5:7, 8:9, drop=FALSE]
  b <- as.matrix(dx[5:7, 8:9])
})

collect()



### --------------------------------------
module("Extraction: 2 columns")

x <- matrix(1:100, ncol=2)
dx <- as.ddmatrix(x)

test("[1:3, 1:2]", {
  a <- x[1:3, 1:2, drop=FALSE]
  b <- as.matrix(dx[1:3, 1:2])
})

test("[1, 1]", {
  a <- x[1, 1, drop=FALSE]
  b <- as.matrix(dx[1, 1])
})

test("[1:2, 1:2]", {
  a <- x[1:2, 1:2, drop=FALSE]
  b <- as.matrix(dx[1:2, 1:2])
})

test("[1:10, -1]", {
  a <- x[1:10, -1, drop=FALSE]
  b <- as.matrix(dx[1:10, -1])
})

collect()



### --------------------------------------
module("Extraction: 2 rows")

x <- matrix(1:100, nrow=2)
dx <- as.ddmatrix(x)

test("[1, 1:2]", {
  a <- x[1, 1:2, drop=FALSE]
  b <- as.matrix(dx[1, 1:2])
})

test("[1, 1]", {
  a <- x[1, 1, drop=FALSE]
  b <- as.matrix(dx[1, 1])
})

test("[1:2, 1:2]", {
  a <- x[1:2, 1:2, drop=FALSE]
  b <- as.matrix(dx[1:2, 1:2])
})

test("[-1, 1:10]", {
  a <- x[-1, 1:10, drop=FALSE]
  b <- as.matrix(dx[-1, 1:10])
})

collect()



finalize()

