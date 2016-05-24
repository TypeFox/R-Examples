library(bfp)

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

exhaustive <- BayesMfp (y ~ bfp (x1, max=1) + bfp(x2, max=1),
                        data = covariateData,
                        priorSpecs =
                        list (a = 3.5,
                              modelPrior="flat"),
                        method = "exhaustive",
                        nModels = 100
                        )           
attr(exhaustive, "logNormConst")
summary(exhaustive)

truedata <- as.data.frame(exhaustive)
truedata

post <- exp(truedata$logMargLik + truedata$logPrior)
normConst <- sum(post)
post / normConst

set.seed(2)
simulation <- BayesMfp (y ~ bfp (x1, max=1) + bfp(x2, max=1),
                        data = covariateData,
                        priorSpecs =
                        list (a = 3.5,
                              modelPrior="flat"),
                        method = "sampling",
                        nModels = 100,
                        chainlength=1000000
                        )           
attr(simulation, "logNormConst")
summary(simulation)

p1 <- posteriors(exhaustive)
p2 <- posteriors(simulation, 2)

plot(log(p1), log(p2))
abline(0, 1)

d1 <- as.data.frame(exhaustive)
d2 <- as.data.frame(simulation)

head(d1[, -c(7,8)])
head(d2[, -c(2, 8, 9)])

shouldbeZero <- d1[, -c(7,8)] - d2[, -c(2, 8, 9)]
shouldbeZero <- max(abs(unlist(shouldbeZero)))
stopifnot(all.equal(shouldbeZero, 0))


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

index <- 1L

stopifnot(all.equal(getLogPrior(dependent[index]),
                    dependent[[index]]$logP))
