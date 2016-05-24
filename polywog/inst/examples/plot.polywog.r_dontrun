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

## All univariate relationships
plot(fit1, n = 20)

## Predicted prestige across occupational categories
plot(fit1, which = "type",
     control.plot = list(xlab = "occupational category"))

## Predicted prestige by education across occupational categories
plot(fit1, which = c("education", "type"), n = 20)

## Joint effect of education and income
plot(fit1, which = c("education", "income"), n = 10)

## Bring up interactive menu
\dontrun{
plot(fit1, ask = TRUE)

  # displays menu:
  # Select one or two variable numbers (separated by spaces), or 0 to exit:

  # 1: education
  # 2: income
  # 3: type
}
