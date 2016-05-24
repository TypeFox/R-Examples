###
### pow2.R  +++ Test suite +++
###


test.pow2 <- function(input, expected) {
    output <- do.call(getFromNamespace("pow2", "pracma"), input)
    identical(output, expected)
}

pow2.expected.00 <- 0
pow2.expected.m1 <- -0.5
pow2.expected.f  <- c(1, 2, 4, 8)
pow2.expected.fe <- c(0, 0.5, -8, 24)
pow2.expected.c  <- c(2^(1i), 2^(-1i))

test.pow2(list(f=0, e=0), pow2.expected.00)
test.pow2(list(f=-1, e=-1), pow2.expected.m1)
test.pow2(list(f=c(0,1,2,3)), pow2.expected.f)
test.pow2(list(f=c(0,1,-2,3), e=c(0,-1,2,3)), pow2.expected.fe)
test.pow2(list(f=c(1i, -1i)), pow2.expected.c)
