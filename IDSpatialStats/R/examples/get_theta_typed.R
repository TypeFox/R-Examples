data(DengueSimulationR02)

r.max<-seq(20,1000,20)
r.min<-seq(0,980,20)

#Lets see if there's a difference in spatial dependence by time case occurs
type<-2-(DengueSimR02[,"time"]<120)
tmp<-cbind(DengueSimR02,type=type)

typed.theta.R01<-get.theta.typed(tmp,typeA=2,typeB=2,r=r.max,r.low=r.min)
