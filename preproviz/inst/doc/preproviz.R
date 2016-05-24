## ----fig.width=10, fig.height=8, message = FALSE, warnings = FALSE-------
library(preproviz)
demoiris <- iris
demoiris[10:20,1] <- NA
a <- preproviz(demoiris)
plotBAR(a)

## ----fig.width=10, fig.height=8, message = FALSE, warnings = FALSE-------
plotHEATMAP(a)
plotVARCLUST(a)

## ----fig.width=10, fig.height=8, message = FALSE, warnings = FALSE-------
plotCMDS(a)

## ----fig.width=10, fig.height=8, message = FALSE, warnings = FALSE-------
plotVARIMP(a)

