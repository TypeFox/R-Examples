## ---- echo = FALSE, warning = FALSE, message = FALSE---------------------
require(GoodmanKruskal)
require(MASS)
require(car)

## ---- echo = TRUE--------------------------------------------------------
GKtau(Cars93$Manufacturer, Cars93$Cylinders)

## ---- echo = TRUE--------------------------------------------------------
GKtau(Cars93$Manufacturer, Cars93$Origin)

## ---- echo = TRUE, fig.width = 7.5, fig.height = 6-----------------------
varSet1 <- c("Manufacturer", "Origin", "Cylinders", "EngineSize", "Passengers")
CarFrame1 <- subset(Cars93, select = varSet1)
GKmatrix1 <- GKtauDataframe(CarFrame1)
plot(GKmatrix1)

## ---- echo = TRUE, fig.width = 7.5, fig.height = 6-----------------------
varSet2 <- c("Manufacturer", "Origin", "Cylinders", "EngineSize", "Passengers", "Make")
CarFrame2 <- subset(Cars93, select = varSet2)
GKmatrix2 <- GKtauDataframe(CarFrame2)
plot(GKmatrix2)

## ---- echo = TRUE--------------------------------------------------------
str(Greene)

## ---- echo = TRUE, fig.width = 7.5, fig.height = 6-----------------------
GKmatrix3 <- GKtauDataframe(Greene)
plot(GKmatrix3)

## ---- echo = TRUE--------------------------------------------------------
table(Greene$language, Greene$location)

## ---- echo = TRUE--------------------------------------------------------
table(mtcars$cyl, mtcars$vs)

## ---- echo = TRUE, fig.width = 7.5, fig.height = 6, fig.cap = "Goodman-Kruskal tau matrix for the mtcars dataframe."----
GKmat <- GKtauDataframe(mtcars)
plot(GKmat, diagSize = 0.8)

## ---- echo = TRUE--------------------------------------------------------
groupedMpg <- GroupNumeric(mtcars$mpg, n = 5)
groupedDisp <- GroupNumeric(mtcars$disp, n = 5)
groupedHp <- GroupNumeric(mtcars$hp, n = 5)
groupedDrat <- GroupNumeric(mtcars$drat, n = 5)
groupedWt <- GroupNumeric(mtcars$wt, n = 5)
groupedQsec <- GroupNumeric(mtcars$qsec, n = 5)
groupedMtcars <- mtcars
groupedMtcars$mpg <- NULL
groupedMtcars$groupedMpg <- groupedMpg
groupedMtcars$disp <- NULL
groupedMtcars$groupedDisp <- groupedDisp
groupedMtcars$hp <- NULL
groupedMtcars$groupedHp <- groupedHp
groupedMtcars$drat <- NULL
groupedMtcars$groupedDrat <- groupedDrat
groupedMtcars$wt <- NULL
groupedMtcars$groupedWt <- groupedWt
groupedMtcars$qsec <- NULL
groupedMtcars$groupedQsec <- groupedQsec

## ---- echo = TRUE, fig.width = 7.5, fig.height = 6, fig.cap = "Goodman-Kruskal tau matrix for the mtcars dataframe."----
GKmat2 <- GKtauDataframe(groupedMtcars)
plot(GKmat2, diagSize = 0.8)

## ---- echo = TRUE--------------------------------------------------------
table(groupedQsec, mtcars$vs)

