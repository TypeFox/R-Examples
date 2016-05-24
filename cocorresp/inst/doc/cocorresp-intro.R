## ----setup, echo = FALSE, results = "hide", message = FALSE--------------
knitr::opts_chunk$set(message = FALSE, warning = FALSE)
library("cocorresp")

## ----load-data-----------------------------------------------------------
data(beetles)
## log transform the beetle data
beetles <- log1p(beetles)
data(plants)

## ----fit-symcoca, message = TRUE-----------------------------------------
bp.sym <- coca(beetles ~ ., data = plants, method = "symmetric")

## ----data-process, message = TRUE----------------------------------------
beetles <- beetles[, colSums(beetles) > 0]
plants <- plants[, colSums(plants) > 0]
bp.sym <- coca(beetles ~ ., data = plants, method = "symmetric")

## ----print-symcoca-------------------------------------------------------
bp.sym

## ----screeplot-symcoca---------------------------------------------------
screeplot(bp.sym)

## ----plot-symcoca, fig.width = 14, fig.height = 7, fig.show = "hold"-----
layout(matrix(1:2, ncol = 2))
plot(bp.sym, which = "response", main = "Beetles")
plot(bp.sym, which = "predictor", main = "Plants")
layout(1)

