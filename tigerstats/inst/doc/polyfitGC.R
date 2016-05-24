## ----include=FALSE-------------------------------------------------------
require(mosaicData)
require(knitr)
opts_chunk$set(tidy=FALSE,fig.width=5.5,fig.height=3.5)

## ----message=FALSE-------------------------------------------------------
require(tigerstats)

## ----eval=FALSE----------------------------------------------------------
#  data(mtcars)
#  View(mtcars)
#  help(mtcars)

## ----eval=FALSE----------------------------------------------------------
#  require(mosaicData)
#  data(RailTrail)
#  View(RailTrail)
#  help(RailTrail)

## ------------------------------------------------------------------------
polyfitGC(mpg ~ wt, data = mtcars, degree=2)

## ----eval=FALSE----------------------------------------------------------
#  polyfitGC(mpg ~ wt, data = mtcars, degree=2,graph=FALSE)

## ------------------------------------------------------------------------
WeightEff <- polyfitGC(mpg ~ wt, data = mtcars,degree=2)

## ------------------------------------------------------------------------
predict(WeightEff, x = 3)

## ------------------------------------------------------------------------
predict(WeightEff,x=3,level=0.95)

## ----fig.width=4.5,fig.height=6------------------------------------------
plot(WeightEff)

## ------------------------------------------------------------------------
polyfitGC(mpg ~ wt, data = mtcars, check = TRUE)

## ------------------------------------------------------------------------
lmGC(volume~avgtemp,data=RailTrail,check=TRUE)

## ------------------------------------------------------------------------
polyfitGC(volume~avgtemp,data=RailTrail,check=TRUE)

