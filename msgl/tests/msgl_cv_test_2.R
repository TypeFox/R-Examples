library(msgl)

# warnings = errors
options(warn=2)

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

fit.cv <- msgl.cv(x, classes, alpha = .5, lambda = lambda, sparse.data = TRUE, standardize = FALSE, max.threads = threads)

err.count <- colSums(fit.cv$classes != classes)

if(err.count[1] < 50 | err.count[25] > 30) stop()
