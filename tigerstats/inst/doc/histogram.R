## ----include=FALSE-------------------------------------------------------
require(knitr)
opts_chunk$set(tidy=FALSE,fig.width=4.5,fig.height=3.5)

## ----message=FALSE-------------------------------------------------------
require(tigerstats)

## ----eval=FALSE----------------------------------------------------------
#  require(lattice)

## ----eval=FALSE----------------------------------------------------------
#  require(devtools)
#  install_github(repo="homerhanumat/tigerstats")

## ----message=FALSE-------------------------------------------------------
require(mosaic)

## ----histtutfastest------------------------------------------------------
histogram(~fastest,data=m111survey,
       type="density",
       xlab="speed (mph)",
       main="Fastest Speed Ever Driven")

## ------------------------------------------------------------------------
favstats(~fastest,data=m111survey)

## ------------------------------------------------------------------------

## ----histtutfastestbreaks10----------------------------------------------
histogram(~fastest,data=m111survey,
       type="density",
       xlab="speed (mph)",
       main="Fastest Speed Ever Driven",
       breaks=c(60,70,80,90,100,110,120,130,
                140,150,160,170,180,190))

## ----histtutfastestbreaks10seq-------------------------------------------
histogram(~fastest,data=m111survey,
       type="density",
       xlab="speed (mph)",
       main="Fastest Speed Ever Driven",
       breaks=seq(from=60,to=190,by=10))

## ----histtutfastestsexcond-----------------------------------------------
histogram(~fastest|sex,data=m111survey,
       type="density",
       xlab="speed (mph)",
       main="Fastest Speed Ever Driven,\nby Sex",
       breaks=seq(from=60,to=190,by=10))

## ----histtutfastestsexcond2----------------------------------------------
histogram(~fastest|sex,data=m111survey,
       type="density",
       xlab="speed (mph)",
       main="Fastest Speed Ever Driven,\nby Sex",
       breaks=seq(from=60,to=190,by=10),
       layout=c(1,2))

