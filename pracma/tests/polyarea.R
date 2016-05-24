###
### polyarea.R  +++ Test suite +++
###


test.polyarea <- function(input, expected) {
   output <- do.call(getFromNamespace("polyarea", "pracma"), input)
   identical(output, expected)
}

polyarea.expected.empty <- 0
polyarea.expected.gen1 <- 3.5
polyarea.expected.gen2 <- 4
polyarea.expected.mtrx <- c(4, 4)
polyarea.expected.cmpl <- 0.5

test.polyarea(list(x=c(), y=c()), polyarea.expected.empty)
test.polyarea(list(x=c(0,2,2,1,0), y=c(0,-1,2,1,1)), polyarea.expected.gen1)
test.polyarea(list(x=matrix(c(1,1,3,3,1), 5, 1),
                   y=matrix(c(1,3,3,1,1), 5, 1)),
              polyarea.expected.gen2)
test.polyarea(list(x=matrix(c(1,3,3,1,1,1,3,3), 4, 2),
                   y=matrix(c(1,1,3,3,1,3,3,1), 4, 2)),
              polyarea.expected.mtrx)
test.polyarea(list(x=c(0,1,1,0), y=c(0,0,1i,0)), polyarea.expected.cmpl)
