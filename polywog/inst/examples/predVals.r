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

## Predicted prestige across occupational categories
predVals(fit1, "type")

## Predicted prestige by education
predVals(fit1, "education", n = 10)

## Plotting
pred_income <- predVals(fit1, "income", n = 10)
plot(pred_income)
