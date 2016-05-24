###
### eig.R  +++ Test suite +++
###


test.eig <- function(input, expected) {
   output <- do.call(getFromNamespace("eig", "pracma"), input)
   identical(output, expected)
}

eig.expected.empty <- matrix(0, nrow=0, ncol=0)
eig.expected.singl <- 1
eig.expected.mat1  <- c(2, 0)
eig.expected.mat2  <- c(1+1i, 1-1i)
eig.expected.mat3  <- c(1, -1)

test.eig(list(a=c()), eig.expected.empty)
test.eig(list(a=c(1)), eig.expected.singl)
test.eig(list(a=matrix(c(1,-1,-1,1), 2, 2)), eig.expected.mat1)
test.eig(list(a=matrix(c(1,1,-1,1), 2, 2)), eig.expected.mat2)
test.eig(list(a=matrix(c(0,1i,-1i,0), 2, 2)), eig.expected.mat3)
