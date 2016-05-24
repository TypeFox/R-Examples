#########################
## Stability selection ##
#########################

library(hdi)

set.seed(123)

x <- matrix(rnorm(100*1000), nrow = 100, ncol = 1000)
y <- x[,1] + x[,2] + rnorm(100)

fit.stab <- stability(x, y, EV = 1)
fit.tmp <- stability(x, y, EV = 1, verbose = TRUE) ## to check verbose

fit.stab
fit.stab$freq[1:10] ## selection frequency of the first 10 predictors

ncores <- if(.Platform$OS.type == "windows") 1 else getOption("mc.cores", 2L)
fit.stab2 <- stability(x, y, EV = 1, parallel = TRUE, ncores = ncores,
                       verbose = TRUE)
fit.stab2
fit.stab2$freq[1:10]
 
stopifnot(all.equal(fit.stab$select, fit.stab2$select))







