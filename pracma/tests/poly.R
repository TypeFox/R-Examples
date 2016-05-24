###
### Poly.R  +++ Test suite +++
###


test.Poly <- function(input, expected) {
   output <- do.call(getFromNamespace("Poly", "pracma"), input)
   identical(output, expected)
}

Poly.expected.empty  <- 1
Poly.expected.1 <- c(1, -6, 11, -6)
#Poly.expected.2 <- error
Poly.expected.3 <- c(1, 0, 0, 0, -1)
Poly.expected.4 <- c(1, -10, 35, -50, 24)
Poly.expected.5 <- c(1, -4, 6, -4, 1)
Poly.expected.6 <- c(1, -5)

test.Poly(list(x=c()), Poly.expected.empty)
test.Poly(list(x=c(1,2,3)), Poly.expected.1)
#test.Poly(list(x=matrix(1:6, 2, 3)), Poly.expected.2)
test.Poly(list(x=c(1,-1,1i,-1i)), Poly.expected.3)
test.Poly(list(x=c(1,2,3,4)), Poly.expected.4)
test.Poly(list(x=diag(4)), Poly.expected.5)
test.Poly(list(x=5), Poly.expected.6)
