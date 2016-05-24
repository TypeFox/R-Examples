library(testthat)
library(Causata)
library(glmnet)
library(XML)

equals <- testthat::equals

context("ToPmml")

# testing data
# creat a data frame with two independent variables and a dependent variable
set.seed(65432)
np <- 10000
# x1 is continuous with some missing values
x1.temp <- rnorm(np)
x1 <- x1.temp
x1[x1>2] <- NA
# x2 is a factor with some missing values
# start with an integer vector that will be used to compute y
x2.int <- sample(-5:5, np, replace=TRUE, prob=c(0.02,0.02,0.02, rep(0.94/8, 8)))
x2.temp <- x2.int
x2.temp[x2.temp == -5] <- NA
x2 <- as.factor(x2.temp)

# binary dependent variable y is a function of x1, x2, and noise
y.temp <- abs(x1.temp)*10 + x2.int + rnorm(np, sd=2)
y <- rep(0,np)
y[y.temp > 15] <- 1

# preprocess data
# using a Discretize transformation
varname <- 'x1__AP'
dvname <- 'y__AP'
df <- data.frame(x1__AP=x1, x2__AP=x2, y__AP=y)
causataData <- CausataData(df, dependent.variable=dvname)
causataData <- ReplaceOutliers(causataData, varname, lowerLimit=-3.5, upperLimit=3.5)
# set breakpoints
breaks <- c(-3.5, -2, -1, 0, 1, 2, 3.5)
fiv <- cut(causataData$df$x1__AP, breaks, include.lowest=TRUE)
# Set weight of evidence, the discrete values that continuous values are mapped to
#woe <- Woe(fiv, causataData$df[[dvname]])
#   [-3.5,-2]      (-2,-1]       (-1,0]        (0,1]        (1,2]      (2,3.5]        BLANK 
# 5.347107531  0.003030305 -4.323101956 -4.360874260 -0.013072082  9.210340372  4.787491743 
woe.levels <- c(5.3, 0, -4.3, -4.4, -0.01, 9.2, 4.8)
# Discretize
causataData <- Discretize(causataData, varname, breaks, woe.levels)

# replace missing values
causataData <- CleanNaFromContinuous(causataData)
causataData <- CleanNaFromFactor(causataData)

# create a model
formula.obj <- formula('y__AP ~ x1__AP + x2__AP')
x.matrix <- model.matrix(formula.obj, data=causataData$df)
cv.glmnet.obj <- cv.glmnet(x.matrix, y=df$y__AP, family='binomial', alpha=1, nfolds=3)

# generate string of PMML
model.def <- ModelDefinition(cv.glmnet.obj, causataData,
  formula.obj, cv.glmnet.obj$lambda.1se )
variable.def <- VariableDefinition(
  name = "score-test-model",
  display.name = "Score: Test Model",
  description = "A logistic regression model.",
  author = "Test User" )
pmml.obj <- ToPmml(model.def, variable.def)

# extract the coefficients from PMML
coefs_pmml <- xpathSApply(pmml.obj, "//x:NumericPredictor", xmlGetAttr, name='coefficient', namespaces='x')
# convert string coefs to numeric
coefs_pmml <- as.numeric(coefs_pmml)
# extract coef names
names(coefs_pmml) <- xpathSApply(pmml.obj, "//x:NumericPredictor", xmlGetAttr, name='name',        namespaces='x')
# replace $ with __ in coef names
names(coefs_pmml) <- gsub('\\$', '__', names(coefs_pmml))
# we have all coefficients except the intercept, add that
coefs_pmml <- c(coefs_pmml, as.numeric(xpathApply(pmml.obj, "//x:RegressionTable", xmlGetAttr, name='intercept', namespaces='x')[1]))
names(coefs_pmml)[length(coefs_pmml)] <- '(Intercept)'

# extract coefficients from glmnet model
coefmatrix <- coef(cv.glmnet.obj)
coefs_glmnet <- as.numeric(coefmatrix)
idx.nonzero <- abs(coefs_glmnet) > 0
coefs_glmnet <- coefs_glmnet[idx.nonzero] # get only nonzero coefficients
names(coefs_glmnet) <- rownames(coefmatrix)[idx.nonzero]

#
# run tests
#

test_that("PMML contains sections for Header, DataDictionary, RegressionModel", 
  expect_that(
    sort(names(xmlChildren(pmml.obj))),
    equals( sort(c('Header', 'DataDictionary', 'RegressionModel')) )
  )
)

test_that("PMML RegressionModel section contains MiningSchema, LocalTransformations, RegressionTable",
  expect_that(
    sort(names(xmlChildren(xmlChildren(pmml.obj)$RegressionModel))),
    equals( sort(c("MiningSchema", "LocalTransformations", "RegressionTable")))
  )
)

test_that("Coefficients in PMML match those in glmnet",
  expect_that(
    # coefficients are ordered by name before the comparison
    coefs_glmnet[order(names(coefs_glmnet))], equals(coefs_pmml[order(names(coefs_pmml))])
  )
)
