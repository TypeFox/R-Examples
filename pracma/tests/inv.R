###
### eig.R  +++ Test suite +++
###


test.inv <- function(input, expected) {
   output <- do.call(getFromNamespace("inv", "pracma"), input)
   identical(output, expected)
}

inv.expected.empty <- matrix(0, nrow=0, ncol=0)
inv.expected.singl <- matrix(Inf, 2, 2)
inv.expected.mat1  <- matrix(c(3,-3,1, -3,5,-2, 1,-2,1), 3, 3)

test.inv(list(a=c()), inv.expected.empty)
test.inv(list(a=matrix(1, 2, 2)), inv.expected.singl)
test.inv(list(a=matrix(c(1,1,1, 1,2,3, 1,3,6), 3, 3)), inv.expected.mat1)
