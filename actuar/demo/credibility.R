### ===== actuar: An R Package for Actuarial Science =====
###
### Demo of the credibility theory facilities provided by actuar
###
### AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>

require(actuar)


## The package provides the famous data set of Hachemeister (1975) as
## a matrix of 5 lines (one for each state) and 25 columns (the state
## number, 12 periods of ratios, 12 periods of corresponding weights).
data(hachemeister)
hachemeister

## Fitting of a Buhlmann model to the Hachemeister data set using
## function 'cm'. The interface of the function is similar to 'lm'.
fit <- cm(~state, hachemeister, ratios = ratio.1:ratio.12)
fit                          # print method
summary(fit)                 # more information
fit$means                    # (weighted) averages
fit$weights                  # total weights
fit$unbiased                 # unbiased variance estimators
predict(fit)                 # credibility premiums

## Fitting of a Buhlmann-Straub model require weights. Here, iterative
## estimators of the variance components are used.
fit <- cm(~state, hachemeister, ratios = ratio.1:ratio.12,
          weights = weight.1:weight.12, method = "iterative")
summary(fit)
predict(fit)

## Simulation of a three level hierarchical portfolio.
nodes <- list(sector = 2, unit = c(3, 4),
              contract = c(10, 5, 8, 5, 7, 11, 4), year = 6)
mf <- expression(sector = rexp(2),
                 unit = rgamma(sector, 0.1),
                 contract = rgamma(unit, 1),
                 year = rpois(weights * contract))
ms <- expression(sector = rnorm(2, sqrt(0.1)),
                 unit = rnorm(sector, 1),
                 contract = NULL,
                 year = rlnorm(unit, 1))
wijkt <- runif(50, 2, 10)
wijkt <- runif(300, rep(0.5 * wijkt, each = 6), rep(1.5 * wijkt, each = 6))
pf <- simul(nodes, model.freq = mf, model.sev = ms, weights = wijkt)

## Fitting of a hierarchical model to the portfolio simulated above.
DB <- cbind(weights(pf, prefix = "weight."),
            aggregate(pf, classif = FALSE) / weights(pf, classif = FALSE))

fit <- cm(~sector + sector:unit + sector:unit:contract,
          data = DB, ratios = year.1:year.6,
          weights = weight.year.1:weight.year.6)
fit
predict(fit)                          # credibility premiums
predict(fit, levels = "unit")         # unit credibility premiums only
summary(fit)                          # portfolio summary
summary(fit, levels = "unit")         # unit portfolio summary only

## Fitting of Hachemeister regression model with intercept at time origin.
fit <- cm(~state, hachemeister, ratios = ratio.1:ratio.12,
          weights = weight.1:weight.12, regformula = ~time,
          regdata = data.frame(time = 1:12))
summary(fit, newdata = data.frame(time = 13)) # 'newdata' is the future value of regressor
predict(fit, newdata = data.frame(time = 13))

## Position the intercept at the barycenter of time.
fit <- cm(~state, hachemeister, ratios = ratio.1:ratio.12,
          weights = weight.1:weight.12, regformula = ~time,
          regdata = data.frame(time = 1:12), adj.intercept = TRUE)
summary(fit, newdata = data.frame(time = 13))
