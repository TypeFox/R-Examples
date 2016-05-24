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
#  data(fuel)
#  View(fuel)
#  help(fuel)

## ----eval=FALSE----------------------------------------------------------
#  require(mosaicData)
#  data(RailTrail)
#  View(RailTrail)
#  help(RailTrail)

## ------------------------------------------------------------------------
lmGC(mpg ~ wt, data = mtcars)

## ------------------------------------------------------------------------
WeightEff <- lmGC(mpg ~ wt, data = mtcars)

## ------------------------------------------------------------------------
predict(WeightEff, x = 3)

## ------------------------------------------------------------------------
predict(WeightEff,x=3,level=0.95)

## ------------------------------------------------------------------------
lmGC(mpg ~ wt, data = mtcars, graph = TRUE)

## ----fig.width=4.5,fig.height=6------------------------------------------
plot(WeightEff)

## ----fig.width=4.5,fig.height=6------------------------------------------
effSpeed <- lmGC(efficiency ~ speed, data = fuel)
plot(effSpeed)

## ------------------------------------------------------------------------
xyplot(efficiency~speed,data=fuel,
        xlab="speed (kilometers/hour",
        ylab="fuel effiency (liters/100km",
        pch=19,col="blue",type=c("p","r"),
        main="Speed and Fuel Efficiency\nfor a British Ford Escort")

## ------------------------------------------------------------------------
lmGC(mpg ~ wt, data = mtcars, check = TRUE)

## ------------------------------------------------------------------------
lmGC(volume~avgtemp,data=RailTrail,check=TRUE)

