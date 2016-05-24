library(msgl)

### Basic tests

data(SimData)
x <- sim.data$x
classes <- sim.data$classes

set.seed(100L)

lambda <- msgl.lambda.seq(x, classes, alpha = .5, d = 25L, lambda.min = 0.05, standardize = TRUE)

if(sgl.c.config()$omp.supported) {
	threads = 2L
} else {
	threads = 1L
}

# This should show a warning
fit.cv <- msgl.cv(x, classes, alpha = .5, fold = 11L, lambda = lambda, standardize = TRUE, max.threads = threads)

err.count <- colSums(fit.cv$classes != classes)

if(err.count[1] < 80 | err.count[25] > 45) stop()
