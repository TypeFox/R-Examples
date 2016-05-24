## Using occupational prestige data
data(Prestige, package = "car")
Prestige <- transform(Prestige, income = income / 1000)

## Fit a polywog model with bootstrap iterations
## (note: using low convergence threshold to shorten computation time of the
## example, *not* recommended in practice!)
set.seed(22)
fit1 <- polywog(prestige ~ education + income + type,
                data = Prestige,
                degree = 2,
                boot = 5,
                thresh = 1e-4)

## Basic information
print(fit1)
summary(fit1)

## See how fitted values change with education holding all else fixed
predVals(fit1, "education", n = 10)

## Plot univariate relationships
plot(fit1)

## Use SCAD instead of adaptive LASSO
fit2 <- update(fit1, method = "scad", thresh = 1e-3)
cbind(coef(fit1), coef(fit2))
