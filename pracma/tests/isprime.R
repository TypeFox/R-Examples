###
### ISPRIME.R  +++ Test suite +++
###


test.isprime <- function(input, expected) {
    output <- do.call(getFromNamespace("isprime", "pracma"), input)
    identical(output, expected)
}

isprime.expected.n1 <- 0
isprime.expected.n2 <- 1
isprime.expected.n100  <- 
    matrix(c(0, 1, 1, 0, 1, 0, 1, 0, 0, 0,
             1, 0, 1, 0, 0, 0, 1, 0, 1, 0,
             0, 0, 1, 0, 0, 0, 0, 0, 1, 0,
             1, 0, 0, 0, 0, 0, 1, 0, 0, 0,
             1, 0, 1, 0, 0, 0, 1, 0, 0, 0,
             0, 0, 1, 0, 0, 0, 0, 0, 1, 0,
             1, 0, 0, 0, 0, 0, 1, 0, 0, 0,
             1, 0, 1, 0, 0, 0, 0, 0, 1, 0,
             0, 0, 1, 0, 0, 0, 0, 0, 1, 0,
             0, 0, 0, 0, 0, 0, 1, 0, 0, 0 ),
           nrow=10, ncol=10, byrow=TRUE)

test.isprime(list(x=1), isprime.expected.n1)
test.isprime(list(x=2), isprime.expected.n2)
test.isprime(list(x=matrix(1:100, 10, 10, byrow=TRUE)),
    isprime.expected.n100)
