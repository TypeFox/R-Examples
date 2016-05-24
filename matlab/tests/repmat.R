###
### $Id: repmat.R 51 2014-02-05 21:22:28Z plroebuck $
###


##-----------------------------------------------------------------------------
test.repmat <- function(input, expected) {
    output <- do.call(getFromNamespace("repmat", "matlab"), input)
    identical(output, expected)
}

X.scalar <- 1
repmat.expected.ones.3x3 <- matlab::ones(3)
repmat.expected.ones.4x2 <- matlab::ones(4, 2)
test.repmat(list(A = X.scalar, n = 3), repmat.expected.ones.3x3)
test.repmat(list(A = X.scalar, m = c(4, 2)), repmat.expected.ones.4x2)
test.repmat(list(A = X.scalar, m = 4, n = 2), repmat.expected.ones.4x2)
test.repmat(list(A = X.scalar, n = matlab::size(repmat.expected.ones.4x2)),
            repmat.expected.ones.4x2)

X.mat <- matlab::eye(2)
repmat.expected.pat1 <- matrix(c(rep(c(1,0), times=3), rep(c(0,1), times=3)),
                               nrow = 6,
                               ncol = 6)
test.repmat(list(A = X.mat, m = 3), repmat.expected.pat1)
test.repmat(list(A = X.mat[, 1], c(1, 1)), X.mat[, 1])

X.vec <- as.numeric(1:8)
repmat.expected.pat2 <- matrix(rep(X.vec, 5), nrow = 5, byrow = TRUE)
test.repmat(list(A = as.matrix(X.vec), m = c(5, 1)), repmat.expected.pat2)

X.str <- "value"
repmat.expected.str <- matrix(rep(X.str, 4), nrow = 1)
test.repmat(list(A = X.str, m = 1, n = 4), repmat.expected.str)

