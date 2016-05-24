library(msgl)

# warnings = errors
options(warn=2)

### Basic tests

data(SimData)
x <- sim.data$x
classes <- sim.data$classes

lambda <- msgl.lambda.seq(x, classes, alpha = .5, d = 25L, lambda.min = 0.05, standardize = TRUE)

test <- replicate(2, 1:20, simplify = FALSE)
train <- lapply(test, function(s) (1:length(classes))[-s])

if(sgl.c.config()$omp.supported) {
	threads = 2L
} else {
	threads = 1L
}

fit.sub <- msgl.subsampling(x, classes, alpha = .5, lambda = lambda, training = train, test = test, max.threads = threads)
if(!all(fit.sub$classes[[1]] == fit.sub$classes[[2]])) stop()
if(min(Err(fit.sub, type="count")) > 15) stop()
