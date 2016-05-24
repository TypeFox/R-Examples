## ---- message = FALSE, results = "hide"----------------------------------
# library(devtools)
# install_github("jslefche/piecewiseSEM")
library(piecewiseSEM)

## ------------------------------------------------------------------------
data(shipley2009)

## ---- message = FALSE, results = "hide"----------------------------------
# Load required libraries for linear mixed effects models
library(lme4)
library(nlme)

# Load example data from package
data(shipley2009)

# Create list of models corresponding to SEM
shipley2009.modlist = list(

  lme(DD ~ lat, random = ~ 1 | site / tree, na.action = na.omit, 
  data = shipley2009),
  
  lme(Date ~ DD, random = ~ 1 | site / tree, na.action = na.omit, 
  data = shipley2009),
  
  lme(Growth ~ Date, random = ~ 1 | site / tree, na.action = na.omit, 
  data = shipley2009),
  
  glmer(Live ~ Growth + (1 | site) + (1 | tree), 
  family = binomial(link = "logit"), data = shipley2009) 
  
  )

## ------------------------------------------------------------------------
sem.fit(shipley2009.modlist, shipley2009, .progressBar = FALSE)

## ------------------------------------------------------------------------
sem.model.fits(shipley2009.modlist)

## ------------------------------------------------------------------------
(coef.table = sem.coefs(shipley2009.modlist, shipley2009))

## ------------------------------------------------------------------------
prod(coef.table$estimate)

## ------------------------------------------------------------------------
sem.coefs(shipley2009.modlist, shipley2009, standardize = "scale")

## ------------------------------------------------------------------------
sem.plot(shipley2009.modlist, shipley2009, standardize = "scale")

## ---- fig.width = 4, fig.height = 4--------------------------------------
# Create new data for predictions
shipley2009.new = data.frame(
  DD = seq(min(shipley2009$DD, na.rm = TRUE), 
           max(shipley2009$DD, na.rm = TRUE), 
           by = 0.01)
)

# Generate predictions
shipley2009.new.pred = sem.predict(shipley2009.modlist, shipley2009.new)
head(shipley2009.new.pred)

# Plot predicted fit
with(shipley2009, plot(Date ~ DD))
lines(shipley2009.new.pred$Date.fit ~ shipley2009.new.pred$DD, lwd = 2, col = "red")

# Generate predictions with standard errors (based on fixed effects only)
shipley2009.new.pred = sem.predict(shipley2009.modlist, shipley2009.new, sefit = TRUE)

# Add 95% confidence bands (roughly 2 * SE)
lines(shipley2009.new.pred$DD, 
      shipley2009.new.pred$Date.fit + 2 * shipley2009.new.pred$Date.se.fit, 
      lwd = 1.5, lty = 2, col = "red")

lines(shipley2009.new.pred$DD, 
      shipley2009.new.pred$Date.fit - 2 * shipley2009.new.pred$Date.se.fit,
      lwd = 1.5, lty = 2, col = "red")

## ------------------------------------------------------------------------
(lavaan.model = sem.lavaan(shipley2009.modlist, shipley2009))


## ---- fig.width = 4, fig.height = 4--------------------------------------
# Create fake data
dat = data.frame(
  y = runif(100),
  x2 = runif(100),
  x3 = runif(100)
)

dat$x1 = dat$y + runif(100, 0, 0.5)

# Create model
model = lm(y ~ x1 + x2 + x3, dat)

# Look at effect of X1 on Y given X2 and X3
partial.resid(y ~ x1, model, dat, return.data.frame = FALSE)

