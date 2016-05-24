data(nhl)

fit<-survfit(Surv(time,status)~1,data=nhl)

plotbands(fit)
