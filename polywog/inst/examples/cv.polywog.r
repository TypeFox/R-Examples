## Using occupational prestige data
data(Prestige, package = "car")
Prestige <- transform(Prestige, income = income / 1000)

## Examine degrees 1 through 4
## (note: using low convergence threshold to shorten computation time of the
## example, *not* recommended in practice!)
set.seed(39)
cv1 <- cv.polywog(prestige ~ education + income + type,
                  data = Prestige,
                  degrees.cv = 1:4,
                  nfolds = 10,
                  thresh = 1e-4)
print(cv1)

## Extract best model and bootstrap
fit1 <- cv1$polywog.fit
fit1 <- bootPolywog(fit1, nboot = 5)
summary(fit1)
