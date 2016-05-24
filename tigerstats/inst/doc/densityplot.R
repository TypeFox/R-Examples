## ----include=FALSE-------------------------------------------------------
require(knitr)
opts_chunk$set(tidy=FALSE,fig.width=5.5,fig.height=3.5)

## ----message=FALSE-------------------------------------------------------
require(tigerstats)

## ----dentutfastest-------------------------------------------------------
densityplot(~fastest,data=m111survey,
       xlab="speed (mph)",
       main="Fastest Speed Ever Driven")

## ----dentutfastestnorug--------------------------------------------------
densityplot(~fastest,data=m111survey,
       xlab="speed (mph)",
       main="Fastest Speed Ever Driven",
       plot.points=FALSE)

## ----dentutfastestsexcond------------------------------------------------
densityplot(~fastest|sex,data=m111survey,
       xlab="speed (mph)",
       main="Fastest Speed Ever Driven,\nby Sex")

## ----dentutfastestsexcond2-----------------------------------------------
densityplot(~fastest|sex,data=m111survey,
       xlab="speed (mph)",
       main="Fastest Speed Ever Driven,\nby Sex",
       layout=c(1,2))

## ----dentutfastestsexgroups----------------------------------------------
densityplot(~fastest,data=m111survey,
       groups=sex,
       xlab="speed (mph)",
       main="Fastest Speed Ever Driven,\nby Sex",
       plot.points=FALSE,
       auto.key=TRUE)

## ----dentutkkard---------------------------------------------------------
densityplot(~kkardashtemp,data=imagpop,
            plot.points=FALSE)

## ----dentutkkardfromto---------------------------------------------------
densityplot(~kkardashtemp,data=imagpop,
            plot.points=FALSE,
            from=0,to=100)

