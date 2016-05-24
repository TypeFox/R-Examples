## ----include=FALSE-------------------------------------------------------
require(knitr)
opts_chunk$set(tidy=FALSE,fig.width=5.5,fig.height=3.5)

## ----message=FALSE-------------------------------------------------------
require(tigerstats)

## ----xytutfastestgpa-----------------------------------------------------
xyplot(fastest~GPA,data=m111survey,
       xlab="grade point average",
       ylab="speed (mph)",
       main="Fastest Speed Ever Driven,\nby Grade Point Average")

## ----xytutfastestgpareg--------------------------------------------------
xyplot(fastest~GPA,data=m111survey,
       xlab="grade point average",
       ylab="speed (mph)",
       main="Fastest Speed Ever Driven,\nby Grade Point Average",
       type=c("p","r"))

## ----xytutfastestgpapoints-----------------------------------------------
xyplot(fastest~GPA,data=m111survey,
       xlab="grade point average",
       ylab="speed (mph)",
       main="Fastest Speed Ever Driven,\nby Grade Point Average",
       pch=19,col="blue")

## ----eval=FALSE----------------------------------------------------------
#  colors()

