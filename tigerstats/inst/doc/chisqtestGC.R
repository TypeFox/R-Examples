## ----include=FALSE-------------------------------------------------------
require(knitr)
require(tigerstats)
opts_chunk$set(tidy=FALSE,fig.width=4.5,fig.height=3.5)

## ------------------------------------------------------------------------
chisqtestGC(~sex+seat,data=m111survey)

## ------------------------------------------------------------------------
SexSeat <- xtabs(~sex+seat,data=m111survey)
SexSeat

## ----eval=FALSE----------------------------------------------------------
#  chisqtestGC(SexSeat)

## ------------------------------------------------------------------------
MySexSeat <- rbind(female=c(19,16,5),male=c(8,16,7))
colnames(MySexSeat) <- c("front","middle","back")

## ------------------------------------------------------------------------
MySexSeat

## ----eval=FALSE----------------------------------------------------------
#  chisqtestGC(MySexSeat)

## ----echo=FALSE----------------------------------------------------------
set.seed(12345)

## ------------------------------------------------------------------------
chisqtestGC(~weather+crowd.behavior,data=ledgejump,
             simulate.p.value="fixed",B=2500)

## ------------------------------------------------------------------------
chisqtestGC(~sex+seat,data=m111survey,
             simulate.p.value="random",B=2500)

## ------------------------------------------------------------------------
chisqtestGC(~sex+seat,data=m111survey,
             simulate.p.value=TRUE,B=2500)

## ----chisqtutnosim, fig.width=4,fig.height=3.5,fig.cap="Graph of P-value, no simulation"----
chisqtestGC(~sex+seat,data=m111survey,graph=TRUE)

## ----chisqtutsim,fig.width=4,fig.height=3.5,fig.cap="Graph of P-value, with simulation"----
chisqtestGC(~sex+seat,data=m111survey,
             simulate.p.value="random",B=2500,graph=TRUE)

## ----chisqtutseatgraph---------------------------------------------------
chisqtestGC(~seat,data=m111survey,
            p=c(1/3,1/3,1/3),
            graph=TRUE)

## ------------------------------------------------------------------------
Seat <- c(front=27,middle=32,back=12)
Seat

## ------------------------------------------------------------------------
chisqtestGC(Seat,p=c(1/3,1/3,1/3))

## ----chisqtutseatsim-----------------------------------------------------
set.seed(678910)
chisqtestGC(~seat,data=m111survey,
            p=c(1/3,1/3,1/3),
            simulate.p.value=TRUE,
            B=2500,graph=TRUE)

## ------------------------------------------------------------------------
chisqtestGC(~sex+seat,data=m111survey,
            verbose=FALSE)

