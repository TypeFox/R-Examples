data(DengueSimulationR01)

r.max<-seq(20,1000,20)
r.min<-seq(0,980,20)

#Lets see if there's a difference in spatial dependence by time case occurs
type<-2-(DengueSimR01[,"time"]<75)
tmp<-cbind(DengueSimR01,type=type)

typed.theta.R01<-get.theta.typed(tmp,typeA=1,typeB=2,r=r.max,r.low=r.min)
typed.theta.type.null<-get.theta.typed.permute(tmp,typeA=1,typeB=2,r=r.max,r.low=r.min,permutations=100)
