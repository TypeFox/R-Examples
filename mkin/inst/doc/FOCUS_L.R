## ---- include = FALSE----------------------------------------------------
library(knitr)
opts_chunk$set(tidy = FALSE, cache = TRUE)

## ------------------------------------------------------------------------
library("mkin")
FOCUS_2006_L1 = data.frame(
  t = rep(c(0, 1, 2, 3, 5, 7, 14, 21, 30), each = 2),
  parent = c(88.3, 91.4, 85.6, 84.5, 78.9, 77.6, 
             72.0, 71.9, 50.3, 59.4, 47.0, 45.1,
             27.7, 27.3, 10.0, 10.4, 2.9, 4.0))
FOCUS_2006_L1_mkin <- mkin_wide_to_long(FOCUS_2006_L1)

## ------------------------------------------------------------------------
m.L1.SFO <- mkinfit("SFO", FOCUS_2006_L1_mkin, quiet=TRUE)
summary(m.L1.SFO)

## ----fig.width=7, fig.height = 5-----------------------------------------
plot(m.L1.SFO)

## ----fig.width=7, fig.height = 5-----------------------------------------
mkinresplot(m.L1.SFO, ylab = "Observed", xlab = "Time")

## ------------------------------------------------------------------------
m.L1.FOMC <- mkinfit("FOMC", FOCUS_2006_L1_mkin, quiet=TRUE)
summary(m.L1.FOMC, data = FALSE)

## ------------------------------------------------------------------------
FOCUS_2006_L2 = data.frame(
  t = rep(c(0, 1, 3, 7, 14, 28), each = 2),
  parent = c(96.1, 91.8, 41.4, 38.7,
             19.3, 22.3, 4.6, 4.6,
             2.6, 1.2, 0.3, 0.6))
FOCUS_2006_L2_mkin <- mkin_wide_to_long(FOCUS_2006_L2)

## ------------------------------------------------------------------------
m.L2.SFO <- mkinfit("SFO", FOCUS_2006_L2_mkin, quiet=TRUE)
summary(m.L2.SFO)

## ----fig.height = 8------------------------------------------------------
par(mfrow = c(2, 1))
plot(m.L2.SFO)
mkinresplot(m.L2.SFO)

## ----fig.height = 8------------------------------------------------------
m.L2.FOMC <- mkinfit("FOMC", FOCUS_2006_L2_mkin, quiet = TRUE)
par(mfrow = c(2, 1))
plot(m.L2.FOMC)
mkinresplot(m.L2.FOMC)
summary(m.L2.FOMC, data = FALSE)

## ----fig.height = 5------------------------------------------------------
m.L2.DFOP <- mkinfit("DFOP", FOCUS_2006_L2_mkin, quiet = TRUE)
plot(m.L2.DFOP)

## ----fig.height = 5------------------------------------------------------
m.L2.DFOP <- mkinfit("DFOP", FOCUS_2006_L2_mkin, 
  parms.ini = c(k1 = 1, k2 = 0.01, g = 0.8),
  quiet=TRUE)
plot(m.L2.DFOP)
summary(m.L2.DFOP, data = FALSE)

## ------------------------------------------------------------------------
FOCUS_2006_L3 = data.frame(
  t = c(0, 3, 7, 14, 30, 60, 91, 120),
  parent = c(97.8, 60, 51, 43, 35, 22, 15, 12))
FOCUS_2006_L3_mkin <- mkin_wide_to_long(FOCUS_2006_L3)

## ----fig.height = 5------------------------------------------------------
m.L3.SFO <- mkinfit("SFO", FOCUS_2006_L3_mkin, quiet = TRUE)
plot(m.L3.SFO)
summary(m.L3.SFO)

## ----fig.height = 5------------------------------------------------------
m.L3.FOMC <- mkinfit("FOMC", FOCUS_2006_L3_mkin, quiet = TRUE)
plot(m.L3.FOMC)
summary(m.L3.FOMC, data = FALSE)

## ----fig.height = 5------------------------------------------------------
m.L3.DFOP <- mkinfit("DFOP", FOCUS_2006_L3_mkin, quiet = TRUE)
plot(m.L3.DFOP)
summary(m.L3.DFOP, data = FALSE)

## ------------------------------------------------------------------------
FOCUS_2006_L4 = data.frame(
  t = c(0, 3, 7, 14, 30, 60, 91, 120),
  parent = c(96.6, 96.3, 94.3, 88.8, 74.9, 59.9, 53.5, 49.0))
FOCUS_2006_L4_mkin <- mkin_wide_to_long(FOCUS_2006_L4)

## ----fig.height = 5------------------------------------------------------
m.L4.SFO <- mkinfit("SFO", FOCUS_2006_L4_mkin, quiet = TRUE)
plot(m.L4.SFO)
summary(m.L4.SFO, data = FALSE)

## ----fig.height = 5------------------------------------------------------
m.L4.FOMC <- mkinfit("FOMC", FOCUS_2006_L4_mkin, quiet = TRUE)
plot(m.L4.FOMC)
summary(m.L4.FOMC, data = FALSE)

