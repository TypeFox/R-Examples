library(parfm)
data(asthma)
head(asthma)
asthma <- asthma[asthma$Fevent==0,]

set.seed(1)
data <- asthma[sample(1:nrow(asthma), 200),]

parfm(Surv(Begin, End, Status)~Drug, cluster="Patid", 
      data=data, dist="weibull", frailty="gamma")

data$time <- data$End - data$Begin

parfm(Surv(time, Status)~Drug, cluster="Patid", data=data,
      dist="weibull", frailty="gamma")
