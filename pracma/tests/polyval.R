###
### polyval.R  +++ Test suite +++
###


test.polyval <- function(input, expected) {
   output <- do.call(getFromNamespace("polyval", "pracma"), input)
   identical(output, expected)
}

polyval.expected.empty1 <- c()
polyval.expected.empty2 <- c(0, 0)
polyval.expected.vec <- c(3, 1, 1, 3, 7)
polyval.expected.mat <- matrix(c(1, 4, 9, 16), nrow=2, ncol=2)

test.polyval(list(p=c(1,1), x=c()), polyval.expected.empty1)
test.polyval(list(p=c(), x=c(1,1)), polyval.expected.empty2)
test.polyval(list(p=c(1,1,1), x=-2:2), polyval.expected.vec)
test.polyval(list(p=c(1,0,0), x=matrix(1:4, 2, 2)), polyval.expected.mat)
