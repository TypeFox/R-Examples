###
### $Id: isprime.R 51 2014-02-05 21:22:28Z plroebuck $
###


##-----------------------------------------------------------------------------
test.isprime <- function(input, expected) {
    output <- do.call(getFromNamespace("isprime", "matlab"), input)
    identical(output, expected)
}

isprime.expected.n1 <- 0
isprime.expected.n2 <- 1
isprime.expected.n100 <- c(0, 1, 1, 0, 1, 0, 1, 0, 0, 0,
                           1, 0, 1, 0, 0, 0, 1, 0, 1, 0,
                           0, 0, 1, 0, 0, 0, 0, 0, 1, 0,
                           1, 0, 0, 0, 0, 0, 1, 0, 0, 0,
                           1, 0, 1, 0, 0, 0, 1, 0, 0, 0,
                           0, 0, 1, 0, 0, 0, 0, 0, 1, 0,
                           1, 0, 0, 0, 0, 0, 1, 0, 0, 0,
                           1, 0, 1, 0, 0, 0, 0, 0, 1, 0,
                           0, 0, 1, 0, 0, 0, 0, 0, 1, 0,
                           0, 0, 0, 0, 0, 0, 1, 0, 0, 0)
isprime.expected.10x10 <- matrix(isprime.expected.n100,
                                 nrow=10, ncol=10, byrow=TRUE)

test.isprime(list(x=1), isprime.expected.n1)
test.isprime(list(x=2), isprime.expected.n2)
test.isprime(list(x=1:100), isprime.expected.n100)
test.isprime(list(x=matrix(1:100, nrow=10, 10, byrow=TRUE)),
                            isprime.expected.n100)

