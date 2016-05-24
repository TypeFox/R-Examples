###
### $Id: sum.R 51 2014-02-05 21:22:28Z plroebuck $
###


##-----------------------------------------------------------------------------
test.sum <- function(input, expected) {
    output <- do.call(getFromNamespace("sum", "matlab"), input)
    identical(all.equal(output,
                        expected,
                        tolerance = 0.0001),
              TRUE)
}

X.vec <- 1:9
sum.expected.vec <- 45

cat("vector test", "\n")
test.sum(list(x = X.vec, na.rm = FALSE), sum.expected.vec)

X.mat <- matrix(X.vec, 3, 3, byrow = TRUE)
sum.expected.mat.by.col <- c(12, 15, 18)
sum.expected.mat.by.row <- c(6, 15, 24)

cat("matrix test", "\n")
test.sum(list(x = X.mat, na.rm = FALSE), sum.expected.mat.by.col)
test.sum(list(x = t(X.mat), na.rm = FALSE), sum.expected.mat.by.row)

X.log <- c(TRUE, TRUE, FALSE, TRUE)
sum.expected.log <- 3

cat("logical test", "\n")
test.sum(list(x = X.log, na.rm = FALSE), sum.expected.log)

