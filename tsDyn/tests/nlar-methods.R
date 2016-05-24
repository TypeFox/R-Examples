library(tsDyn)

x <- log10(lynx)

### Estimate models
mod <- list()
mod[["linear"]] <- linear(x, m=2)
mod[["setar"]] <- setar(x, m=2, thDelay=1)
mod[["lstar"]] <- lstar(x, m=2, thDelay=1)
mod[["aar"]] <- aar(x, m=2)


### Extract methods
sapply(mod, AIC)
sapply(mod, BIC)
sapply(mod, mse)
sapply(mod, MAPE)

sapply(mod, coef)
sapply(mod, function(x) head(residuals(x)))


lapply(mod, predict, n.ahead=10)


### Pred Roll, acc_stat:
x_small <- x[1:100]
mod_small <- list()
mod_small[["linear"]] <- linear(x_small, m=2)
mod_small[["setar"]] <- setar(x_small, m=2, thDelay=1, th=getTh(mod[["setar"]]))
mod_small[["lstar"]] <- lstar(x_small, m=2, thDelay=1, th=getTh(mod[["lstar"]]), gamma=coef(mod[["lstar"]])["gamma"])
mod_small[["aar"]] <- aar(x_small, m=2)

pred_rolls_1 <- lapply(mod_small, predict_rolling, n.ahead=1, newdata=x[101:114])
sapply(pred_rolls_1, function(x) x$pred[[1]])
sapply(pred_rolls_1, accuracy_stat)


pred_rolls_12 <- lapply(mod_small, predict_rolling, n.ahead=1:2, newdata=x[101:114])
sapply(pred_rolls_12, function(x) x$pred[[1]])
lapply(pred_rolls_12, accuracy_stat)

