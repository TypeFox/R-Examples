library(condvis)
data(powerplant)
## fit a model
models <- list(
    lm = lm(PE ~ ., data = powerplant),
    lmquad = lm(PE ~ . + I(AT^2), data = powerplant))
## visualise sections along 'AT'
ceplot(data = powerplant, model = models, S = "AT")