###
### NTHROOT.R  +++ Test suite +++
###


test.nthroot <- function(input, expected) {
    output <- do.call(getFromNamespace("nthroot", "pracma"), input)
    identical(output, expected)
}

nthroot.expected.n0 <- c(-1.2)
nthroot.expected.n1 <- c(1, -2, 3)
nthroot.expected.n2 <- c(1,  0, 3)
nthroot.expected.n3 <- c(1, -2, 3)

test.nthroot(list(x=-1.2^5, n=5), nthroot.expected.n0)
test.nthroot(list(x=c(1,-2, 3), n=1), nthroot.expected.n1)
test.nthroot(list(x=c(1, 0, 9), n=2), nthroot.expected.n2)
test.nthroot(list(x=c(1,-8,27), n=3), nthroot.expected.n3)
