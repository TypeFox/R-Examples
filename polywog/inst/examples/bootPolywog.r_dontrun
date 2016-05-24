## Using occupational prestige data
data(Prestige, package = "car")
Prestige <- transform(Prestige, income = income / 1000)

## Fit a polywog model without bootstrapping
## (note: using low convergence threshold to shorten computation time of the
## example, *not* recommended in practice!)
fit1 <- polywog(prestige ~ education + income + type,
                data = Prestige,
                degree = 2,
                thresh = 1e-4)
summary(fit1)

## Bootstrap the fitted model
fit2 <- bootPolywog(fit1, nboot = 5)
summary(fit2)

## Example of parallel processing on Mac/Unix via 'doMC'
\dontrun{
library(doMC)
registerDoMC()

fit2 <- bootPolywog(fit1, nboot = 100, .parallel = TRUE)
}

## Example of parallel processing on Windows via 'doSMP'
\dontrun{
library(doSMP)
w <- startWorkers()
registerDoSMP(w)

fit2 <- bootPolywog(fit1, nboot = 100, .parallel = TRUE)

stopWorkers(w)
}
