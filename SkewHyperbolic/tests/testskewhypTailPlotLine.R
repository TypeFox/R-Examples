require(SkewHyperbolic)
## source("../R/skewhypTailPlotLine.R")
### Draw tail plot of some data
param <- c(0,1,2,10)
set.seed(101, kind = "default", normal.kind = "default")
x <- rskewhyp(500, param = param)
tailPlot(x)
### Add skew hyperbolic t-distribution line
skewhypTailPlotLine(x, param = param)
### Parameters from fit may look better
paramFit <- skewhypFit(x, plots = FALSE)$param
tailPlot(x)
skewhypTailPlotLine(x, param = param)
skewhypTailPlotLine(x, param = paramFit, col = "steelblue")

### Left tail example
tailPlot(x, side = "l")
### Add skew hyperbolic t-distribution line
skewhypTailPlotLine(x, param = paramFit, side = "l")
### Log scale on both axes
tailPlot(x, side = "r", log = "xy")
### Add skew hyperbolic t-distribution line
skewhypTailPlotLine(x, param = paramFit, side = "r")

