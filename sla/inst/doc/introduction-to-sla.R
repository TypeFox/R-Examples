## ----set_opts, echo = FALSE----------------------------------------------
options(show.signif.stars = FALSE)

## ------------------------------------------------------------------------
library(sla)

## ------------------------------------------------------------------------
data(eqslo)
eqslo

## ------------------------------------------------------------------------
obj <- sla(eqslo)

## ------------------------------------------------------------------------
obj

## ------------------------------------------------------------------------
summary(obj)

## ----abcd_plot, fig.cap = "Four models of interest fit to the eqslo data: (A) Individual intercepts and individual slopes, (B) Single intercept and single slope, (C) Individual intercepts and single slope, (D) Single intercept and individual slopes."----
par(mfrow = c(2, 2))
plot(obj, modelType2Plot = 'A')
plot(obj, modelType2Plot = 'B')
plot(obj, modelType2Plot = 'C')
plot(obj, modelType2Plot = 'D')

## ------------------------------------------------------------------------
library(sla) 
library(MASS)
data(whiteside)

## ------------------------------------------------------------------------
wsobj <- sla(whiteside)

## ------------------------------------------------------------------------
wsobj

## ------------------------------------------------------------------------
attributes(wsobj)

## ------------------------------------------------------------------------
summary(wsobj$Mod.A)

## ------------------------------------------------------------------------
summary(wsobj)

## ------------------------------------------------------------------------
anova(wsobj$Mod.C, wsobj$Mod.A)

## ----a_plot, fig.cap = "Model A: Individual intercepts and individual slopes for the whiteside data."----
plot(wsobj, mod = 'A')

