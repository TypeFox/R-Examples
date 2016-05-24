###
### $Id: flipud.R 51 2014-02-05 21:22:28Z plroebuck $
###


##-----------------------------------------------------------------------------
test.flipud <- function(input, expected) {
    output <- do.call(getFromNamespace("flipud", "matlab"), input)
    identical(output, expected)
}

X.mat <- matrix(1:6, 3, 2)
flipud.expected.X.mat <- matrix(c(3:1, 6:4), 3, 2)

test.flipud(list(X.mat), flipud.expected.X.mat)

X.vec <- seq(1, 9, by = 2)
flipud.expected.X.vec <- seq(9, 1, by = -2)

test.flipud(list(X.vec), flipud.expected.X.vec)

## Since R does not really have concept of row/column vectors
## identical(test.fliplr(X.vec, flipud.expected.X.vec),
##           test.flipud(X.vec, flipud.expected.X.vec))

