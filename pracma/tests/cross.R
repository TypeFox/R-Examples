###
### cross.R  +++ Test suite +++
###


test.cross <- function(input, expected) {
   output <- do.call(getFromNamespace("cross", "pracma"), input)
   identical(output, expected)
}

cross.expected.1 <- c(-3, 6, -3)
#cross.expected. <- 
#cross.expected. <- 
#cross.expected. <- 

test.cross(list(x=c(1, 2, 3), y=c(4, 5, 6)), cross.expected.1)
#test.cross(list(), cross.expected.)
#test.cross(list(), cross.expected.)
#test.cross(list(), cross.expected.)
