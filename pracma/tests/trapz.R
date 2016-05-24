###
### trapz.R  +++ Test suite +++
###


test.trapz <- function(input, expected) {
   output <- do.call(getFromNamespace("trapz", "pracma"), input)
   identical(output, expected)
}

trapz.expected.empty1 <- 0
trapz.expected.empty2 <- 0
trapz.expected.gen1 <- 12
trapz.expected.gen2 <- 6
trapz.expected.cmpl1 <- 0+0.5i
trapz.expected.cmpl2 <- 0+0.5i

test.trapz(list(x=c()), trapz.expected.empty1)
test.trapz(list(x=c(), y=c()), trapz.expected.empty2)
test.trapz(list(x=1:5), trapz.expected.gen1)
test.trapz(list(x=seq(0,2,by=0.5), y=1:5), trapz.expected.gen2)
test.trapz(list(x=c(0,1), y=c(0,1i)), trapz.expected.cmpl1)
test.trapz(list(x=c(0,1i), y=c(0,1)), trapz.expected.cmpl2)
