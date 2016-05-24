###
### hilb.R  +++ Test suite +++
###


test.hilb <- function(input, expected) {
   output <- do.call(getFromNamespace("hilb", "pracma"), input)
   identical(output, expected)
}

hilb.expected.m1 <- matrix(NA, nrow=0, ncol=0)
hilb.expected.0  <- matrix(0, nrow=0, ncol=0)
hilb.expected.1  <- matrix(1, nrow=1, ncol=1)
hilb.expected.5  <- 1 / matrix(c(1:5,2:6,3:7,4:8,5:9), nrow=5,ncol=5)

test.hilb(list(n=-1), hilb.expected.m1)
test.hilb(list(n=0), hilb.expected.0)
test.hilb(list(n=1), hilb.expected.1)
test.hilb(list(n=5), hilb.expected.5)
