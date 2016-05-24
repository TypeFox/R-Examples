suppressMessages(library(cobs))
options(digits = 6, warn = 2) ## << all warnings to errors!

## When 'R CMD check'ing, we may want to see exact package information:
sessionInfo() # plus the details of the major dependent packages:
packageDescription("SparseM")
packageDescription("quantreg")
packageDescription("cobs")
##

set.seed(101)
x <- seq(-2,2,length = 100)
y <- x^2+0.5*rnorm(100)
## Constraints -- choosing ones that are true for  f(x) = x^2
PW  <- rbind(
             c(0, -3,9), # f(-3) = 9
             c(0,  3,9), # f(3 ) = 9
             c(2,  0,0)) # f'(0) = 0

mod <- cobs (x,y,constraint = "convex", pointwise = PW)
mod

stopifnot(all.equal(predict(mod, c(-3, 3))[,"fit"], c(9,9), tol = 1e-12))

## derivative 0 at 0 -- we miss a 'deriv = 1' argument [-> see ../TODO]
eps <- 1e-6
stopifnot(abs(diff(predict(mod, c(-eps, eps))[,"fit"])/(2*eps)) < .001 * eps)

