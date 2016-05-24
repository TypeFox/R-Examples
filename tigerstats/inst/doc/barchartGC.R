## ----include=FALSE-------------------------------------------------------
require(knitr)
opts_chunk$set(tidy=FALSE,fig.width=5.5,fig.height=3.5)

## ----message=FALSE-------------------------------------------------------
require(tigerstats)

## ----eval=FALSE----------------------------------------------------------
#  data(m111survey)
#  View(m111survey)
#  help(m111survey)

## ----barchartseatfreq,tidy=FALSE-----------------------------------------
barchartGC(~seat,data=m111survey,type="frequency",
           main="Barchart of Seating Preference",
           xlab="Seating Preference")

## ----barchartseatperc,tidy=FALSE-----------------------------------------
barchartGC(~seat,data=m111survey,type="percent",
           main="Barchart of Seating Preference",
          xlab="Seating Preference")

## ------------------------------------------------------------------------
Seat <- xtabs(~seat,data=m111survey)
Seat

## ----barchartseattab, tidy=FALSE-----------------------------------------
barchartGC(Seat,type="percent",
           main="Barchart of Seating Preference",
           xlab="Seating Preference")

## ----barchartsexseat, tidy=FALSE-----------------------------------------
barchartGC(~sex+seat,data=m111survey,
           type="percent",
           main="Sex and Seating Preference\nat Georgetown",
           xlab="Sex of student",
           sub="Bar color shows seating preference")

## ------------------------------------------------------------------------
SexSeat <- xtabs(~sex+seat,data=m111survey)
SexSeat

## ----barchartsexseattab, tidy=FALSE--------------------------------------
barchartGC(SexSeat,type="percent",
           main="Sex and Seating Preference\nat Georgetown",
           xlab="Sex of student",
           sub="Bar color shows seating preference")

## ------------------------------------------------------------------------
rowPerc(SexSeat)

## ------------------------------------------------------------------------
barchartGC(SexSeat,type="percent",
           main="Sex and Seating Preference\nat Georgetown",
           ylab="Sex of student",
           sub="Bar color shows seating preference",
           flat=TRUE)

## ------------------------------------------------------------------------
barchartGC(SexSeat,type="freq",
           main="Sex and Seating Preference\nat Georgetown",
           xlab="Sex of student",
           sub="Bar color shows seating preference",
           stack=TRUE)

## ------------------------------------------------------------------------
barchartGC(SexSeat,type="freq",
           main="Sex and Seating Preference\nat Georgetown",
           ylab="Sex of student",
           sub="Bar color shows seating preference",
           stack=FALSE,
           horizontal=TRUE)

## ----barchartfastest,tidy=FALSE,fig.width=7,fig.height=5-----------------
barchartGC(~fastest,data=m111survey,
           main="Fastest Speed Ever Driven")

