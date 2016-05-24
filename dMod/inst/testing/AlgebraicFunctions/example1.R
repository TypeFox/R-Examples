# Framework
library(cOde)
library(dMod)

# Plotting
library(ggplot2)
library(ggthemes)

# Optimization
library(trust)


## Skript -------------------------------

# Beobachtungsfunktion ("time" ist reserviert als unabhängige Variable)
y <- Y(c(x = "a*time + b*time^2"), f = NULL)

# Prädiktor
predictor <- matrix(0:10, ncol = 1, dimnames = list(NULL, "time"))
attr(predictor, "sensitivities") <- predictor

# Daten
data <- data.frame(time = 0:10, name = "x", value  = 2*(0:10) + 0.5 * (0:10)^2 + rnorm(11, 0, 3), sigma = 3)
plotData(list(condition1 = data))

# Objective function
obj <- function(p) wrss(res(data, y(predictor, p)))

# Fitten
myfit <- trust(obj, parinit = c(a = 1, b = 1), rinit = 1, rmax = 10)
plotCombined(list(condition1 = y(predictor, myfit$argument)), list(condition1 = data))

# Profile
profiles.approx <- sapply(1:2, function(i) profile.trust(obj, myfit$argument, whichPar = i, limits = c(-10, 10)))
plotProfile(profiles.approx)
plotPaths(profiles.approx)




