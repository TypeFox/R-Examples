###
### $Id: size.R 51 2014-02-05 21:22:28Z plroebuck $
###


##-----------------------------------------------------------------------------
test.size <- function(input, expected) {
    output <- do.call(getFromNamespace("size", "matlab"), input)
    identical(as.integer(output), expected)
}

X.vec <- 2:9
size.expected.X.vec <- c(1, length(X.vec))

test.size(list(X = X.vec), size.expected.X.vec)
test.size(list(X = X.vec, 1), as.integer(1))
test.size(list(X = X.vec, 2), length(X.vec))


X.mat <- matrix(X.vec, 4, 2)
size.expected.X.mat <- dim(X.mat)

test.size(list(X = X.mat), size.expected.X.mat)
test.size(list(X = X.mat, 1), size.expected.X.mat[1])
test.size(list(X = X.mat, 2), size.expected.X.mat[2])
test.size(list(X = X.mat, 3), as.integer(1))	# singleton dimension

X.arr <- array(2:25, c(4, 3, 2))
size.expected.X.arr <- dim(X.arr)

test.size(list(X = X.arr), size.expected.X.arr)
test.size(list(X = X.arr, 1), size.expected.X.arr[1])
test.size(list(X = X.arr, 2), size.expected.X.arr[2])
test.size(list(X = X.arr, 3), size.expected.X.arr[3])

