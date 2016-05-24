###
### $Id: hilb.R 51 2014-02-05 21:22:28Z plroebuck $
###


##-----------------------------------------------------------------------------
test.hilb <- function(input, expected) {
    output <- do.call(getFromNamespace("hilb", "matlab"), input)
    identical(output, expected)
}

hilb.expected.na.nodim   <- matrix(NA, nrow=0, ncol=0)
hilb.expected.zero.nodim <- matrix(0, nrow=0, ncol=0)
hilb.expected.one <- matrix(1, nrow=1, ncol=1)
hilb.expected.five <- 1 / matrix(c(1:5,2:6,3:7,4:8,5:9), nrow=5, ncol=5)

test.hilb(list(n=-1), hilb.expected.na.nodim)
test.hilb(list(n=0), hilb.expected.zero.nodim)
test.hilb(list(n=1), hilb.expected.one)
test.hilb(list(n=5), hilb.expected.five)

