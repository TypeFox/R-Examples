library(bfp)

set.seed(19)

x1 <- rnorm(n=15)
x2 <- rbinom(n=15, size=20, prob=0.5) 
x3 <- rexp(n=15)

y <- rt(n=15, df=2)

## run an exhaustive model space evaluation with a flat model prior and
## a uniform prior (a = 4) on the shrinkage factor t = g/(1 + g):
test <- BayesMfp(y ~ bfp (x2, max = 4) + uc (x1 + x3), nModels = 100,
                 method="exhaustive")

test <- BayesMfp(y ~ bfp(x2, max=1) + uc (x3), nModels = 100,
                 method="exhaustive",
                 priorSpecs=list(a=4, modelPrior="dependent"))
summary(test)
logPriors <- as.data.frame(test)$logPrior

sum(exp(logPriors[c(1)]))
sum(exp(logPriors[c(2, 4, 12:18)]))
sum(exp(logPriors[c(3, 5:11)]))

sum(exp(logPriors[3:4])) / sum(exp(logPriors[3:18]))

test <- BayesMfp(y ~ bfp (x2, max = 4) + uc (x1 + x3), nModels = 100,
                 method="exhaustive", priorSpecs=list(a=4, modelPrior="sparse"))

test <- BayesMfp(y ~ bfp (x2, max = 4) + uc (x1 + x3), nModels = 100,
                 method="exhaustive", priorSpecs=list(a=4, modelPrior="dependent"))

for(index in seq_along(test))
{
    stopifnot(all.equal(getLogPrior(test[index]),
                        test[[index]]$logP))
}


summary(test)
## setting

beta0 <- 1
alpha1 <- 1
alpha2 <- 3
delta1 <- 1

sigma <- 2                              # sigma2 = 4
n <- 15
k <- 2L

## simulate data

set.seed (123)

x <- matrix (runif (n * k, 1, 4), nrow = n, ncol = k) # predictor values
w <- matrix (rbinom (n * 1, size = 1, prob = 0.5), nrow = n, ncol = 1)

x1tr <- alpha1 * x[,1]^2
x2tr <- alpha2 * (x[,2])^(1/2)
w1tr <- delta1 * w[,1]

predictorTerms <-
    x1tr +
    x2tr +
    w1tr

trueModel <- list (powers = list (x1 = 2, x2 = 0.5),
                   ucTerms = as.integer (1)
                   )

covariateData <- data.frame (x1 = x[,1],
                             x2 = x[,2],
                             w = w)

covariateData$y <- predictorTerms + rnorm (n, 0, sigma)
covariateData

## also test the dependent model prior
dependent <- BayesMfp (y ~ bfp (x1, max=1) + bfp(x2, max=1),
                        data = covariateData,
                        priorSpecs =
                        list (a = 3.5,
                              modelPrior="dependent"),
                        method = "exhaustive",
                        nModels = 10000)           
attr(dependent, "logNormConst")

depSum <- as.data.frame(dependent)
depSum

sum(exp(depSum$logPrior))

for(index in seq_along(dependent))
{
    stopifnot(all.equal(getLogPrior(dependent[index]),
                        dependent[[index]]$logP))
}


## and with sampling:
## the frequencies must converge to
## the (normalized) posterior probabilities! 
set.seed(93)
dependent2 <- BayesMfp(y ~ bfp (x1, max=1) + bfp(x2, max=1),
                        data = covariateData,
                        priorSpecs =
                        list (a = 3.5,
                              modelPrior="sparse"),
                        method = "sampling",
                        nModels = 10000L,
                       chainlength=1000000L)      
depSum2 <- as.data.frame(dependent2)
depSum2

