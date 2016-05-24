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
alpha <- 0
d <- 20L
lambda.min <- 0.5
algorithm.config <- sgl.standard.config 

# create data
data <- create.sgldata(x, y, weights, sampleGrouping)
lambda <- sgl_lambda_sequence("sgl_test_dense", "sglOptim", data, covariateGrouping, groupWeights, parameterWeights, alpha = alpha, d = d, lambda.min, algorithm.config)
fit1a.cv <- sgl_cv("sgl_test_dense", "sglOptim", data, covariateGrouping, groupWeights, parameterWeights, alpha, lambda, fold = 2L, cv.indices = list(), max.threads = 1L, algorithm.config)

data <- create.sgldata(x, y, weights, sampleGrouping, sparseX = TRUE)

#Seed used for cv splitting
set.seed(100)

fit1b.cv <- sgl_cv("sgl_test_sparse", "sglOptim", data, covariateGrouping, groupWeights, parameterWeights, alpha, lambda, fold = 2L, cv.indices = list(), max.threads = 1L, algorithm.config)
