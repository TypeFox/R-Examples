## ----include=FALSE-------------------------------------------------------
require(knitr)
opts_chunk$set(tidy=FALSE,fig.width=4.5,fig.height=3.5)

## ----message=FALSE-------------------------------------------------------
require(tigerstats)

## ----bwtutfastest--------------------------------------------------------
bwplot(~fastest,data=m111survey,
       xlab="speed (mph)",
       main="Fastest Speed Ever Driven")

## ----bwtutfastestsex-----------------------------------------------------
bwplot(fastest~sex,data=m111survey,
       ylab="speed (mph)",
       main="Fastest Speed Ever Driven,\nby Sex of Subject")

## ----bwtutfastestsexhoriz------------------------------------------------
bwplot(sex~fastest,data=m111survey,
       xlab="speed (mph)",
       main="Fastest Speed Ever Driven,\nby Sex of Subject")

## ----bwtutbad------------------------------------------------------------
bwplot(fastest~GPA,data=m111survey,
       ylab="speed (mph)",
       xbal="grade-point average",
       main="Fastes Speed Ever Driven,\nby Grade-Point Average")

