###
### compan.R  +++ Test suite +++
###


test.compan <- function(input, expected) {
   output <- do.call(getFromNamespace("compan", "pracma"), input)
   identical(output, expected)
}

compan.expected.empty  <- c()
compan.expected.sngl1  <- c()
compan.expected.sngl2  <- c()
compan.expected.bspl1  <- matrix(c(0, 7, -6, 1, 0, 0, 0, 1, 0),
                                 nrow=3, ncol=3, byrow=TRUE)

test.compan(list(p=c()), compan.expected.empty)
test.compan(list(p=c(0)), compan.expected.sngl1)
test.compan(list(p=c(1)), compan.expected.sngl2)
test.compan(list(p=c(1,0,-7,6)), compan.expected.bspl1)
