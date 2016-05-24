###
### $Id: rot90.R 51 2014-02-05 21:22:28Z plroebuck $
###


##-----------------------------------------------------------------------------
test.rot90 <- function(input, expected) {
    output <- do.call(getFromNamespace("rot90", "matlab"), input)
    identical(output, expected)
}

X.mat <- matrix(1:9, nrow = 3, byrow = TRUE)
rot90.expected.X.mat  <- matrix(c(3:1, 6:4, 9:7), nrow = 3)
rot180.expected.X.mat <- matrix(c(9:7, 6:4, 3:1), nrow = 3, byrow = TRUE)
rot270.expected.X.mat <- matrix(c(7:9, 4:6, 1:3), nrow = 3)
rot360.expected.X.mat <- X.mat

test.rot90(list(A = X.mat, k = 1), rot90.expected.X.mat)
test.rot90(list(A = X.mat, k = 2), rot180.expected.X.mat)
test.rot90(list(A = X.mat, k = 3), rot270.expected.X.mat)
test.rot90(list(A = X.mat, k = 4), rot360.expected.X.mat)

