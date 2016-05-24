library(grplasso)

set.seed(719)

n <- 100
p <- 10

x <- matrix(runif(n * p, min = -2.5, max = 2.5), nrow = n, ncol = p)
y <- 4 * sin(x[,1]) + x[,2]^2 + rnorm(n)

#####################
##                 ## 
## grplasso object ##
##                 ##
#####################

fit <- grplasso(cbind(1,x), y, index = c(NA, rep(1:5, each = 2)),
                lambda = c(400, 200, 100, 50, 25, 10, 5), standardize = FALSE,
                model = LinReg(), control = grpl.control(trace = 0))

##plot(fit, log = "x")

## See whether sub-setting is working correctly
## Can't work with all.equal here because of attributes
if(any(range(coef(fit)[,3:4] - coef(fit[3:4])) != 0))
  stop("Subsetting not working correctly")

## Check whether errors are produced for non-valid examples
m1 <- inherits(try(fit[1:12], silent = TRUE), "try-error")

## Try to remove everything
m2 <- inherits(try(fit[-(1:8)], silent = TRUE), "try-error")

## Stop if any of m1, ... is *FALSE*
stopifnot(m1, m2)


