library(sglOptim)
# warnings = errors
options(warn=2)

data(TestData)
x <- test.data$x
y <- test.data$y
grp <- test.data$grp

weights <- rep(1/nrow(x), nrow(x))
sampleGrouping <- grp
covariateGrouping <- factor(1:ncol(x))
groupWeights <- c(sqrt(length(levels(sampleGrouping))*table(covariateGrouping)))
parameterWeights <-  matrix(1, nrow = length(levels(sampleGrouping)), ncol = ncol(x))
alpha <- 0.5
d <- 20L
lambda.min <- 0.5

# create data
data <- create.sgldata(x, y, weights, sampleGrouping)
lambda <- sgl_lambda_sequence("sgl_test_dense", "sglOptim", data, covariateGrouping, groupWeights, parameterWeights, alpha = alpha, d = d, lambda.min)

# indices
test <- replicate(2, 1:20, simplify = FALSE)
train <- lapply(test, function(s) (1:nrow(x))[-s])

if(sgl.c.config()$omp.supported) {
	threads = 2L
} else {
	threads = 1L
}
	
fit1a.sub <- sgl_subsampling("sgl_test_dense", "sglOptim", data, covariateGrouping, groupWeights, parameterWeights, alpha, lambda, train, test, max.threads = threads)
