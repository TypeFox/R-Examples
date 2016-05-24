num1to100.gen <- setNumericWithRange("Numeric", min = 1, max = 100)
par.gen <- setRefClass("Graph",
                       properties(list(size = "NumericWithMin1Max100")))
pars <- par.gen$new(size = new("NumericWithMin1Max100", 5))
pars$size #current value is 5
try(pars$size <- 300) # out of range error
pars$size <- 10 #works

## Positive Integer
par.gen <- setRefClass("PI", properties(list(size  = "PositiveInteger"),
                                        list(size = PositiveInteger(2))))
obj <- par.gen$new()
## error
try(obj$size <- -1)
obj$size <- 3
