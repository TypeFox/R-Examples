
suppressMessages(library(RcppDE))

## somewhat pathodological example with nuisance parameter mul
Rastrigin <- function(x) {
    mul * (sum(x+2 - 10 * cos(2*pi*x)) + 20)
}

## create a new environment associated with the function
funenv <- environment(fun=Rastrigin)
assign("mul", 2, envir=funenv)        ## set value

out <- DEoptim(Rastrigin, -25, 25,
               control=list(NP=10, trace=0),
               env=funenv)
summary(out)
