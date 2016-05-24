###
### dot.R  +++ Test suite +++
###


test.dot <- function(input, expected) {
   output <- do.call(getFromNamespace("dot", "pracma"), input)
   identical(output, expected)
}

dot.expected.empty <- 0
dot.expected.55 <- 55
dot.expected.t55 <- 55
dot.expected.mm <- c(26, 44)
dot.expected.00 <- 0
dot.expected.neg <- -2

test.dot(list(x=c(), y=c()), dot.expected.empty)
test.dot(list(x=1:5, y=1:5), dot.expected.55)
test.dot(list(x=1:5, y=t(t(1:5))), dot.expected.t55)
test.dot(list(x=matrix(c(1,3,2,4), 2, 2), y=matrix(c(5,7,6,8), 2, 2)),
         dot.expected.mm)
test.dot(list(x=c(0, 0), y=c(1, 2)), dot.expected.00)
test.dot(list(x=c(1, 1), y=c(-1, -1)), dot.expected.neg)
