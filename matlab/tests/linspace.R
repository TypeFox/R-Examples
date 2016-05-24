###
### $Id: linspace.R 51 2014-02-05 21:22:28Z plroebuck $
###


##-----------------------------------------------------------------------------
test.linspace <- function(input, expected) {
    output <- do.call(getFromNamespace("linspace", "matlab"), input)
    identical(output, expected)
}

linspace.expected.100 <- as.numeric(1:100)
linspace.expected.25x5 <- c(0, 6.25, 12.50, 18.75, 25)
linspace.expected.len1 <- 25

test.linspace(list(a = 1, b = 100), linspace.expected.100)
test.linspace(list(a = 0, b = 25, n = 5), linspace.expected.25x5)
test.linspace(list(a = 1, b = 25, n = 1), linspace.expected.len1)

## more rigorously this time
test.linspace(list(a = 0, b = 1, n = 2.5), 0:1)     ## HWB 2011/02/03
test.linspace(list(a = 0, b = 1, n = 3.9), seq(0, 1, length=floor(3)))

