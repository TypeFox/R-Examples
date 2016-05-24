## This is a small simulation study.

## Setting:
beta0 <- 1
alpha1 <- 1
alpha2 <- c (1, 1)
delta1 <- 1

sigma <- 2 ## => sigma2 = 4
n <- 40

## simulate data:
set.seed (123)

## variables (continuous x and binary w)
x <- matrix (runif (n * 3, 1, 4), nrow = n, ncol = 3) 
w <- matrix (rbinom (n * 2, size = 1, prob = 0.5), nrow = n, ncol = 2)

## transformed continuous covariates
x1tr <- alpha1 * x[,1]^2
x2tr <- cbind (x[,2]^-0.5, x[,2]^-0.5 * log (x[,2])) %*% alpha2

summary (x1tr / x2tr)
mean (x1tr) / mean (x2tr)

w1tr <- delta1 * w[,1]

## so the linear predictor is
predictorTerms <-
    x1tr +
    x2tr +
    w1tr

## define the true model
trueModel <- list (powers=
                   list(x.1=2,
                        x.2=c (-0.5, -0.5),
                        x.3=numeric(0)),
                   ucTerms = as.integer (1))

covariateData <- data.frame (x = x, w = w)


## simulate the response from the linear predictor and error standard deviation
covariateData$y <- predictorTerms + rnorm (n, 0, sigma)

## and see where the true model is ranked
system.time(models <- BayesMfp (y ~ bfp (x.1) + bfp (x.2) + bfp (x.3) + uc (w.1) + uc (w.2),
                                data = covariateData,
                                method = "exhaustive"))
ind <- findModel(trueModel, models)

## get the estimated effects in the true model and compare with the true shapes
estimate <- models[ind]
summary (estimate)

plotx1 <- plotCurveEstimate (estimate, "x.1")
with (plotx1, lines (original, original^2, col = "red"))

plotx2 <- plotCurveEstimate (estimate, "x.2")
with (plotx2, lines (original, original^-0.5 * (1 + log (original)), col = "red"))

## what is the order if equal prior probabilities are used for the models?
modelsSummary <- as.data.frame (models)
modelsSummary$rankML <- as.character (rank (-modelsSummary$logMargLik))
modelsSummary[1:20, ]

library (doBy)
orderBy (~ -logMargLik, data = modelsSummary)[1:100,] 


