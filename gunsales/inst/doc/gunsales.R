## ---- setup, include=FALSE-----------------------------------------------
library(gunsales)

## ---- initialDataFake, eval=FALSE----------------------------------------
#  gunsales <- analysis()

## ---- initialData, fig.width=7, echo=FALSE-------------------------------
if (gunsales:::.goodOS()) gunsales <- analysis() else cat("Unsupported platform -- no plots below.")

## ---- basePlotsFake, eval=FALSE------------------------------------------
#  plot_gunsales(gunsales)

## ---- basePlots, fig.width=7, echo=FALSE---------------------------------
if (gunsales:::.goodOS()) plot_gunsales(gunsales)

## ---- ggPlotsFake, eval=FALSE--------------------------------------------
#  ggplot_gunsales(gunsales)

## ---- ggPlots, fig.width=7, echo=FALSE-----------------------------------
if (gunsales:::.goodOS()) ggplot_gunsales(gunsales)

