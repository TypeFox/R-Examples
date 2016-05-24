library(parfm)
data(insem)
head(insem)

insem$TimeMonths <- insem$Time * 12 / 365.25
  
set.seed(1)
insem <- insem[sample(1:nrow(insem), 150),]
parfm(Surv(TimeMonths, Status)~Heifer, cluster="Herd", data=insem,
             dist="weibull", frailty="gamma")
