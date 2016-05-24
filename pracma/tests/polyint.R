###
### polyint.R  +++ Test suite +++
###


test.polyint <- function(input, expected) {
   output <- do.call(getFromNamespace("polyint", "pracma"), input)
   identical(output, expected)
}

polyint.expected.1 <- c(1, 0)
polyint.expected.2 <- c(1/6, 1/5, 1/4, 1/3, 1/2, 1, 1)

test.polyint(list(p=c(1)), polyint.expected.1)
test.polyint(list(p=c(1,1,1,1,1,1), k=1), polyint.expected.2)
