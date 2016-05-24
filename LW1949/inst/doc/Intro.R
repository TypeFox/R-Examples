## ----load_package--------------------------------------------------------
library(LW1949)

## ----dose-effect_data----------------------------------------------------
conc <- c(0.0625, 0.125, 0.25, 0.5, 1, 2, 3)
numtested <- rep(8, 7)
numaffected <- c(1, 4, 4, 7, 8, 8, 8)
mydat <- dataprep(dose=conc, ntot=numtested, nfx=numaffected)

## ----dataprep()_output---------------------------------------------------
mydat

## ----fit_LW--------------------------------------------------------------
intslope <- fitLWauto(mydat)
fLW <- LWestimate(intslope, mydat)

## ----fitLWauto_output----------------------------------------------------
intslope

## ----LWestimate_output---------------------------------------------------
fLW

## ----estimate_LW---------------------------------------------------------
pctaffected <- c(25, 50, 99.9)
predlinear(pctaffected, fLW)

## ----plot_fits, fig.width=5, fig.height=5--------------------------------
plotDELP(mydat)
predLinesLP(fLW)

plotDE(mydat)
predLines(fLW)

