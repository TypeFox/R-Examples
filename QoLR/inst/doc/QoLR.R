## ------------------------------------------------------------------------
library(QoLR)
data(dataqol)
head(dataqol)

## ------------------------------------------------------------------------
score_dataqol=scoring.QLQC30(dataqol,id="Id", time="time")
head(round(score_dataqol))

## ------------------------------------------------------------------------
info=dataqol[,c("Id","time","date","death","Arm")]
dataqol_final=merge(score_dataqol,info,by=c("Id","time"))

## ------------------------------------------------------------------------
dataqol_final=dataqol_final[,c(1:2,18,3:17,19:20)]
head(round(dataqol_final))

## ------------------------------------------------------------------------
dataqol_final=dataqol_final[order(dataqol_final$time),]
dataqol_final=dataqol_final[order(dataqol_final$Id),]
head(round(dataqol_final))

## ------------------------------------------------------------------------
ttd1=TTD(dataqol_final, score="QL", MCID=5)
head(ttd1)

## ------------------------------------------------------------------------
ttd2=TTD(dataqol_final, score="QL", order=1, MCID=5, no_baseline="event",no_follow="event")
head(ttd2)

## ------------------------------------------------------------------------
ttd3=TTD(dataqol_final, score="QL", MCID=5, death="death")
head(ttd3)

## ------------------------------------------------------------------------
ttd4=TTD(dataqol_final, score="QL", MCID=5, death="death", sensitivity=TRUE)
head(ttd4)

## ------------------------------------------------------------------------
ttd5=TTD(dataqol_final, score="QL", MCID=5, ref.init ="best")
head(ttd5)

## ------------------------------------------------------------------------
ttd6=TTD(dataqol_final, score=c("QL","PF","FA"), order=c(1,1,2), MCID=5)
head(ttd6)

## ------------------------------------------------------------------------
tudd1=TUDD(dataqol_final, score="QL", MCID=5)
head(tudd1)

## ------------------------------------------------------------------------
tudd2=TUDD(dataqol_final, score = "QL", MCID = 5, ref.def = "def3")
head(tudd2)

## ------------------------------------------------------------------------
tudd3=TUDD(dataqol_final, score="PF", MCID=c(5,10), sensitivity=T)
head(round(tudd3,2))

## ----plotTTD, fig.width=5.5, fig.align='center'--------------------------
tudd1=TUDD(dataqol_final, score="QL", MCID=5,ref.init="baseline",ref.def="def1")
ttd_1=merge(tudd1,unique(dataqol_final[,c("Id","Arm")]))
plotTTD(ttd_1$time.5.QL,ttd_1$event.5.QL,ttd_1$Arm,nrisk=T,nevent=F,
group.names=c("Arm 1","Arm 2"), t=seq(0,10,2),info=T,pos.info=c(6,0.8),
xlab="time (months)", ylab="probability (%)")

