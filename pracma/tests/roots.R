###
### roots.R  +++ Test suite +++
###


test.roots <- function(input, expected) {
   output <- do.call(getFromNamespace("roots", "pracma"), input)
   identical(output, expected)
}

roots.expected.empty <- matrix(0, nrow=0, ncol=0)
roots.expected.singl <- matrix(0, nrow=0, ncol=0)
roots.expected.bspl1 <-  c(0, 2, -2, 1, -1)  # Matlab: c(0, -2, -1, 1, 2)
                         c(0, 2, -2, 1, -1)
roots.expected.bspl2 <-  c(0.5, -0.2)
roots.expected.bspl3 <- -c(0, 0, -1, 1)

test.roots(list(p=c()), roots.expected.empty)
test.roots(list(p=c(0)), roots.expected.singl)
#test.roots(list(p=c(1,0,-5,0,4,0)), roots.expected.bspl1)  # zapsmall
test.roots(list(p=c(1,-0.3,-0.1)), roots.expected.bspl2)
test.roots(list(p=c(1,0,-1,0,0)), roots.expected.bspl3)
