## ----include=FALSE-------------------------------------------------------
library(knitr)
opts_chunk$set(
concordance=TRUE
)

## ----echo=FALSE, message=FALSE-------------------------------------------
require(knitr)
opts_chunk$set(cache = TRUE, dev = "pdf")

## ----message=F, echo=TRUE------------------------------------------------
# load MVN package
library(MVN)

## ----message=F, eval=FALSE, echo=TRUE------------------------------------
#  # load Iris data
#  data(iris)

## ----message=FALSE, echo=TRUE--------------------------------------------
# setosa subset of the Iris data
setosa <- iris[1:50, 1:4]

## ----"Mardia test", message=FALSE, echo=TRUE-----------------------------
result <- mardiaTest(setosa, qqplot = FALSE)
result

## ----message=FALSE, echo=FALSE-------------------------------------------
tmp <- mardiaTest(setosa, qqplot = FALSE)

## ----"Henze-Zirkler test", message=FALSE, echo=TRUE----------------------
result <- hzTest(setosa, qqplot = FALSE)
result

## ----"Royston test", message=FALSE, echo=FALSE---------------------------
result <- roystonTest(setosa, qqplot = FALSE)
result

## ----echo=FALSE, message=FALSE, fig.width=5, fig.height=5----------------
par(mar=c(4.2,4.1,3,0.2))
result <- roystonTest(setosa, qqplot = TRUE)

## ----qqUniPlot, eval=FALSE, message=FALSE, echo=TRUE---------------------
#  uniPlot(setosa, type = "qqplot") # creates univariate Q-Q plots
#  uniPlot(setosa, type = "histogram") # creates univariate histograms

## ----include=TRUE, echo=FALSE--------------------------------------------
par(cex.main=1)
uniPlot(setosa, type= "qqplot")

## ----include=TRUE, echo=FALSE--------------------------------------------
par(cex.main=1)
uniPlot(setosa, type= "histogram")

## ----eval=FALSE, message=FALSE, echo=TRUE--------------------------------
#  uniNorm(setosa, type = "SW", desc = TRUE)

## ----SWUnivariate, eval=TRUE, message=FALSE, echo=FALSE------------------
setosa = iris[1:50,-5]
uniNorm(setosa, type = "SW", desc = TRUE)

## ----message=FALSE, echo=FALSE, fig.keep='none'--------------------------
setosa3 <- iris[1:50, 1:3] 
mard <- mardiaTest(setosa3, qqplot=FALSE)
hz <- hzTest(setosa3, qqplot=FALSE)
roys <- roystonTest(setosa3, qqplot=FALSE)

## ----eval=FALSE, fig.keep='none', message=FALSE, echo=TRUE---------------
#  setosa2 <- iris[1:50, 1:2]
#  result <- hzTest(setosa2, qqplot=FALSE)
#  mvnPlot(result, type = "persp", default = TRUE) # perspective plot
#  mvnPlot(result, type = "contour", default = TRUE) # contour plot

## ----message=FALSE, echo=FALSE, fig.width=5, fig.height=5----------------
setosa2 <- iris[1:50, 1:2]
result <- hzTest(setosa2, qqplot = FALSE)
par <- par(mar=c(0.4,0.4,0.4,0.4))
mvnPlot(result, type = "persp", default = TRUE)

## ----echo=FALSE, message=FALSE, fig.width=5, fig.height=5----------------
setosa2 <- iris[1:50, 1:2]
result <- hzTest(setosa2)
mvnPlot(result, type = "contour", default = TRUE)

## ----message=FALSE, echo=FALSE, fig.keep='none'--------------------------
setosa2 <- iris[1:50, 1:2] 
mard <- mardiaTest(setosa2, qqplot = FALSE)
hz <- hzTest(setosa2, qqplot = FALSE)
roys <- roystonTest(setosa2, qqplot = FALSE)

## ----mvoutlier, eval=FALSE, echo=TRUE------------------------------------
#  versicolor <- iris[51:100, 1:3]
#  # Mahalanobis distance
#  result <- mvOutlier(versicolor, qqplot = TRUE, method = "quan")
#  # Adjusted Mahalanobis distance
#  result <- mvOutlier(versicolor, qqplot = TRUE, method = "adj.quan")

## ----echo=FALSE, message=FALSE, fig.width=5, fig.height=5----------------
versicolor <- iris[51:100, 1:3]
result <- mvOutlier(versicolor, qqplot = TRUE, method = "quan")

## ----echo=FALSE, message=FALSE, fig.width=5, fig.height=5----------------
versicolor <- iris[51:100, 1:3]
result <- mvOutlier(versicolor, qqplot = TRUE, method = "adj.quan")

