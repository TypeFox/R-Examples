###
### Test with q=1 (each feature has one parameter)
###

library(sglOptim)

# warnings = errors
options(warn=2)

data(TestData)
x <- test.data$x
y <- test.data$y
grp <- test.data$grp

weights <- rep(1/nrow(x), nrow(x))
sampleGrouping <- factor(rep(1, nrow(x))) #Note
covariateGrouping <- factor(1:ncol(x))
groupWeights <- c(sqrt(length(levels(sampleGrouping))*table(covariateGrouping)))
parameterWeights <-  matrix(1, nrow = length(levels(sampleGrouping)), ncol = ncol(x))
d <- 50L
lambda.min <- 2
algorithm.config <- sgl.standard.config 

# To check dimension do
#data <- create.sgldata(x, y, weights, sampleGrouping)
#args <- prepare.args(data, covariateGrouping, groupWeights, parameterWeights, alpha = 0)
#args$block.dim


# Fit tests dense
# create data
data <- create.sgldata(x, y, weights, sampleGrouping)
# alpha = 0
lambda <- sgl_lambda_sequence("sgl_test_dense", "sglOptim", data, covariateGrouping, groupWeights, parameterWeights, alpha = 0, d = d, lambda.min, algorithm.config)
fit1a <- sgl_fit("sgl_test_dense", "sglOptim", data, covariateGrouping, groupWeights, parameterWeights, alpha = 0, lambda, return = 1:length(lambda), algorithm.config)
# alpha = 0.5
lambda <- sgl_lambda_sequence("sgl_test_dense", "sglOptim", data, covariateGrouping, groupWeights, parameterWeights, alpha = 0.5, d = d, lambda.min, algorithm.config)
fit1a <- sgl_fit("sgl_test_dense", "sglOptim", data, covariateGrouping, groupWeights, parameterWeights, alpha = 0.5, lambda, return = 1:length(lambda), algorithm.config)
# alpha = 1
lambda <- sgl_lambda_sequence("sgl_test_dense", "sglOptim", data, covariateGrouping, groupWeights, parameterWeights, alpha = 1, d = d, lambda.min, algorithm.config)
fit1a <- sgl_fit("sgl_test_dense", "sglOptim", data, covariateGrouping, groupWeights, parameterWeights, alpha = 1, lambda, return = 1:length(lambda), algorithm.config)

# Predict test
res1a <- sgl_predict("sgl_test_dense", "sglOptim", fit1a, data)

# Fit tests sparse
data <- create.sgldata(x, y, weights, sampleGrouping, sparseX = TRUE)
fit1b <- sgl_fit("sgl_test_sparse", "sglOptim", data, covariateGrouping, groupWeights, parameterWeights, alpha = 0, lambda, return = 1:length(lambda), algorithm.config)
fit1b <- sgl_fit("sgl_test_sparse", "sglOptim", data, covariateGrouping, groupWeights, parameterWeights, alpha = 0.5, lambda, return = 1:length(lambda), algorithm.config)
fit1b <- sgl_fit("sgl_test_sparse", "sglOptim", data, covariateGrouping, groupWeights, parameterWeights, alpha = 1, lambda, return = 1:length(lambda), algorithm.config)

# Predict test
res1b <- sgl_predict("sgl_test_sparse", "sglOptim", fit1b, data)


if(max(abs(fit1a$beta[[25]]-fit1b$beta[[25]])) > 1e-5) stop()
if(max(abs(res1a$responses$link[[25]]-res1b$responses$link[[25]])) > 1e-5) stop()
