###
### polyder.R  +++ Test suite +++
###


test.polyder <- function(input, expected) {
   output <- do.call(getFromNamespace("polyder", "pracma"), input)
   identical(output, expected)
}

polyder.expected.0 <- 0
polyder.expected.1 <- 0
polyder.expected.5 <- c(4, 3, 2, 1)
polyder.expected.3 <- c(2, 0)
polyder.expected.2 <- c(12, 36, 42, 18)


test.polyder(list(p=c()), polyder.expected.0)
test.polyder(list(p=c(1)), polyder.expected.1)
test.polyder(list(p=c(1,1,1,1,1)), polyder.expected.5)
test.polyder(list(p=c(1,0,0), q=c(0,0,1)), polyder.expected.3)
test.polyder(list(p=c(3,6,9), q=c(1,2,0)), polyder.expected.2)
