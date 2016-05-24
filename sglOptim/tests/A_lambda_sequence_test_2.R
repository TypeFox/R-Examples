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
d <- 100L
lambda.min <- 0.1
algorithm.config <- sgl.standard.config 

# To check dimension do
#data <- create.sgldata(x, y, weights, sampleGrouping)
#args <- prepare.args(data, covariateGrouping, groupWeights, parameterWeights, alpha = 0)
#args$block.dim

# create data
data <- create.sgldata(x, y, weights, sampleGrouping)
lambda_da <- sgl_lambda_sequence("sgl_test_dense", "sglOptim", data, covariateGrouping, groupWeights, parameterWeights, alpha = 0, d = d, lambda.min, algorithm.config)

data <- create.sgldata(x, y, weights, sampleGrouping, sparseX = TRUE)
lambda_sp <- sgl_lambda_sequence("sgl_test_sparse", "sglOptim", data, covariateGrouping, groupWeights, parameterWeights, alpha = 0, d = d, lambda.min, algorithm.config)

if(max(abs(lambda_sp-lambda_da)) > 1e-5) stop()

data <- create.sgldata(x, y, weights, sampleGrouping)
lambda_da <- sgl_lambda_sequence("sgl_test_dense", "sglOptim", data, covariateGrouping, groupWeights, parameterWeights, alpha = 1, d = d, lambda.min, algorithm.config)

data <- create.sgldata(x, y, weights, sampleGrouping, sparseX = TRUE)
lambda_sp <- sgl_lambda_sequence("sgl_test_sparse", "sglOptim", data, covariateGrouping, groupWeights, parameterWeights, alpha = 1, d = d, lambda.min, algorithm.config)

if(max(abs(lambda_sp-lambda_da))  > 1e-5) stop()

data <- create.sgldata(x, y, weights, sampleGrouping)
lambda_da <- sgl_lambda_sequence("sgl_test_dense", "sglOptim", data, covariateGrouping, groupWeights, parameterWeights, alpha = .5, d = d, lambda.min, algorithm.config)

data <- create.sgldata(x, y, weights, sampleGrouping, sparseX = TRUE)
lambda_sp <- sgl_lambda_sequence("sgl_test_sparse", "sglOptim", data, covariateGrouping, groupWeights, parameterWeights, alpha = .5, d = d, lambda.min, algorithm.config)

if(max(abs(lambda_sp-lambda_da))  > 1e-5) stop()
