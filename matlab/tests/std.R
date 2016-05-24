###
### $Id: std.R 51 2014-02-05 21:22:28Z plroebuck $
###


##-----------------------------------------------------------------------------
test.std <- function(input, expected) {
    output <- do.call(getFromNamespace("std", "matlab"), input)
    identical(all.equal(output,
                        expected,
                        tolerance = 0.0001),
              TRUE)
}

X.mat <- matrix(c(1, 5, 9, 7, 15, 22), 2, 3, byrow = TRUE)
std.expected.by.col <- c(4.2426, 7.0711, 9.1924)
std.expected.by.row <- c(4.000, 7.5056)

test.std(list(x = X.mat), std.expected.by.col)
test.std(list(x = t(X.mat)), std.expected.by.row)

