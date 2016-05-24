###
### $Id: rem.R 51 2014-02-05 21:22:28Z plroebuck $
###


##-----------------------------------------------------------------------------
test.rem <- function(input, expected) {
    output <- do.call(getFromNamespace("rem", "matlab"), input)
    identical(output, expected)
}

X.mat <- matrix(1:9, nrow = 3, byrow = TRUE)
rem.expected.X.mat <- matrix(c(3:1, 6:4, 9:7), nrow = 3, byrow = TRUE)
rem.expected.X.mat.Y0 <- matrix(rep(NaN, length(X.mat)), nrow = nrow(X.mat))
rem.expected.X.mat.Y1 <- matrix(rep(1, length(X.mat)), nrow = nrow(X.mat))
rem.expected.X.mat.Y2 <- matrix(rep(c(1, 0), 5)[1:length(X.mat)],
                                nrow = nrow(X.mat))
rem.expected.X.mat.Y3 <- matrix(rep(c(1, 2, 0), nrow(X.mat)),
                                nrow = nrow(X.mat),
                                byrow = TRUE)

test.rem(list(x = X.mat, y = 0), rem.expected.X.mat.Y0)
test.rem(list(x = X.mat, y = 1), rem.expected.X.mat.Y1)
test.rem(list(x = X.mat, y = 2), rem.expected.X.mat.Y2)
test.rem(list(x = X.mat, y = 3), rem.expected.X.mat)


## rem & mod give same results with X, Y having same sign
test.rem(list(x = 5, y = 3), matlab::mod(5, 3))
test.rem(list(x = -5, y = -3), matlab::mod(-5, -3))

## alternate formula used when X, Y having different signs
test.rem(list(x = 5, y = -3), (matlab::mod(5, -3) - -3))
test.rem(list(x = -5, y = 3), (matlab::mod(-5, 3) - 3))

