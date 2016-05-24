###
### $Id: find.R 51 2014-02-05 21:22:28Z plroebuck $
###


##-----------------------------------------------------------------------------
test.find <- function(input, expected) {
    output <- do.call(getFromNamespace("find", "matlab"), input)
    identical(output, expected)
}

X <- c(3, 2, 1, 1, 2, 3)
find.expected.eq.one <- c(3, 4)

test.find(list(X == 1), find.expected.eq.one)

X <- c(1, 0, 4, -3, 0, 0, 0, 8, 6)
find.expected.nonzero <- c(1, 3, 4, 8, 9)
find.expected.gt.two <- c(3, 8, 9)

test.find(list(X), find.expected.nonzero)
test.find(list(X > 2), find.expected.gt.two)

