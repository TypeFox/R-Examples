##### several checks for function elementary_symmetric_functions()
##### 2012/09/07, BAEK


##### preliminaries
library("psychotools")
set.seed(1)


##################################################
##### dichotomous items (zero and first order derivatives)
##################################################
di10 <- runif(10, -2.5, 2.5)
di20 <- runif(20, -2.5, 2.5) 
elementary_symmetric_functions(di10)
elementary_symmetric_functions(di20)

## with 10 items, everything is fine 
all.equal(elementary_symmetric_functions(di10, order = 1, engine = "R", diff = FALSE), elementary_symmetric_functions(di10, order = 1, engine = "C", diff = FALSE)) # R-bin sum vs. C sum
all.equal(elementary_symmetric_functions(di10, order = 1, engine = "R", diff = FALSE), elementary_symmetric_functions(di10, order = 1, engine = "C", diff = TRUE))  # R-bin sum vs. C diff
all.equal(elementary_symmetric_functions(di10, order = 1, engine = "R", diff = TRUE), elementary_symmetric_functions(di10, order = 1, engine = "C", diff = FALSE))  # R-bin diff vs. C sum
all.equal(elementary_symmetric_functions(di10, order = 1, engine = "R", diff = TRUE), elementary_symmetric_functions(di10, order = 1, engine = "C", diff = TRUE))   # R-bin diff vs. C diff
all.equal(elementary_symmetric_functions(as.list(di10), order = 1, engine = "R", diff = FALSE), elementary_symmetric_functions(di10, order = 1, engine = "R", diff = FALSE)) # R-poly sum vs. R-bin sum
all.equal(elementary_symmetric_functions(as.list(di10), order = 1, engine = "R", diff = FALSE), elementary_symmetric_functions(di10, order = 1, engine = "R", diff = TRUE)) # R-poly sum vs. R-bin diff
all.equal(elementary_symmetric_functions(as.list(di10), order = 1, engine = "R", diff = FALSE), elementary_symmetric_functions(di10, order = 1, engine = "R", diff = FALSE)) # R-poly diff vs. R-bin sum
all.equal(elementary_symmetric_functions(as.list(di10), order = 1, engine = "R", diff = FALSE), elementary_symmetric_functions(di10, order = 1, engine = "R", diff = TRUE)) # R-poly diff vs. R-bin diff

## with 20 items we get (rounding) errors for difference algorithm :)
all.equal(elementary_symmetric_functions(di20, order = 1, engine = "R", diff = FALSE), elementary_symmetric_functions(di20, order = 1, engine = "C", diff = FALSE)) # R-bin sum vs. C sum
all.equal(elementary_symmetric_functions(di20, order = 1, engine = "R", diff = FALSE), elementary_symmetric_functions(di20, order = 1, engine = "C", diff = TRUE))  # R-bin sum vs. C diff
all.equal(elementary_symmetric_functions(di20, order = 1, engine = "R", diff = TRUE), elementary_symmetric_functions(di20, order = 1, engine = "C", diff = FALSE))  # R-bin diff vs. C sum
all.equal(elementary_symmetric_functions(di20, order = 1, engine = "R", diff = TRUE), elementary_symmetric_functions(di20, order = 1, engine = "C", diff = TRUE))   # R-bin diff vs. C diff
all.equal(elementary_symmetric_functions(as.list(di20), order = 1, engine = "R", diff = FALSE), elementary_symmetric_functions(di20, order = 1, engine = "R", diff = FALSE)) # R-poly sum vs. R-bin sum
all.equal(elementary_symmetric_functions(as.list(di20), order = 1, engine = "R", diff = FALSE), elementary_symmetric_functions(di20, order = 1, engine = "R", diff = TRUE)) # R-poly sum vs. R-bin diff
all.equal(elementary_symmetric_functions(as.list(di20), order = 1, engine = "R", diff = FALSE), elementary_symmetric_functions(di20, order = 1, engine = "R", diff = FALSE)) # R-poly diff vs. R-bin sum
all.equal(elementary_symmetric_functions(as.list(di20), order = 1, engine = "R", diff = FALSE), elementary_symmetric_functions(di20, order = 1, engine = "R", diff = TRUE)) # R-poly diff vs. R-bin diff

## second order derivatives (only implemented for R-bin sum/diff)
all.equal(elementary_symmetric_functions(di10, order = 2, engine = "R", diff = FALSE), elementary_symmetric_functions(di10, order = 2, engine = "R", diff = TRUE)) # R-bin sum vs. R-bin diff
all.equal(elementary_symmetric_functions(di20, order = 2, engine = "R", diff = FALSE), elementary_symmetric_functions(di10, order = 2, engine = "R", diff = TRUE)) # R-bin sum vs. R-bin diff


##################################################
##### polytomous items (zero and first order derivatives)
##################################################
pi10 <- unname(split(runif(10, -2.5, 2.5), factor(rep(1:5, each = 2))))
pi20 <- unname(split(runif(20, -2.5, 2.5), factor(rep(1:10, each = 2))))
elementary_symmetric_functions(pi10)
elementary_symmetric_functions(pi20)

## with 10 items, everything is fine.
all.equal(elementary_symmetric_functions(pi10, order = 1, engine = "R", diff = FALSE), elementary_symmetric_functions(pi10, order = 1, engine = "C", diff = FALSE)) # R-poly sum vs. C sum
all.equal(elementary_symmetric_functions(pi10, order = 1, engine = "R", diff = FALSE), elementary_symmetric_functions(pi10, order = 1, engine = "C", diff = TRUE))  # R-poly sum vs. C diff
all.equal(elementary_symmetric_functions(pi10, order = 1, engine = "R", diff = TRUE), elementary_symmetric_functions(pi10, order = 1, engine = "C", diff = FALSE))  # R-poly diff vs. C sum
all.equal(elementary_symmetric_functions(pi10, order = 1, engine = "R", diff = TRUE), elementary_symmetric_functions(pi10, order = 1, engine = "C", diff = TRUE))   # R-poly diff vs. C diff

## with 20 items we get (rounding) errors for difference algorithm :)
all.equal(elementary_symmetric_functions(pi20, order = 1, engine = "R", diff = FALSE), elementary_symmetric_functions(pi20, order = 1, engine = "C", diff = FALSE)) # R-poly sum vs. C sum
all.equal(elementary_symmetric_functions(pi20, order = 1, engine = "R", diff = FALSE), elementary_symmetric_functions(pi20, order = 1, engine = "C", diff = TRUE))  # R-poly sum vs. C diff
all.equal(elementary_symmetric_functions(pi20, order = 1, engine = "R", diff = TRUE), elementary_symmetric_functions(pi20, order = 1, engine = "C", diff = FALSE))  # R-poly diff vs. C sum
all.equal(elementary_symmetric_functions(pi20, order = 1, engine = "R", diff = TRUE), elementary_symmetric_functions(pi20, order = 1, engine = "C", diff = TRUE))   # R-poly diff vs. C diff


##################################################
##### check 'check'
##################################################
all.equal(elementary_symmetric_functions(di10, order = 1), elementary_symmetric_functions(as.list(di10), order = 1)) # dichotomous items as list / vector
inherits(try(elementary_symmetric_functions(di10, log = FALSE, engine = "C")), "try-error")                          # input on log scale for dichotmous items not possible (engine C)
inherits(try(elementary_symmetric_functions(pi10, log = FALSE, engine = "R")), "try-error")                          # input on log scale for polytomous items not possible (engine R)
inherits(try(elementary_symmetric_functions(pi10, log = FALSE, engine = "C")), "try-error")                          # input on log scale for polytomous items not possible (engine C)
inherits(try(elementary_symmetric_functions(di10, order = 2, engine = "C")), "try-error")                            # second order derivatives for dichotmous items not possible (engine C)
inherits(try(elementary_symmetric_functions(pi10, order = 2, engine = "R")), "try-error")                          # second order derivatives for polytomous items not possible (engine R)
inherits(try(elementary_symmetric_functions(pi10, order = 2, engine = "C")), "try-error")                          # second order derivatives for polytomous items not possible (engine C)
