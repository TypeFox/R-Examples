data(DengueSimulationR01)

r.max<-seq(20,1000,20)
r.min<-seq(0,980,20)

#Lets see if there's a difference in spatial dependence by time case occurs
type<-2-(DengueSimR01[,"time"]<75)
tmp<-cbind(DengueSimR01,type=type)

typed.theta.bs<-get.theta.typed.bootstrap(tmp,typeA=1,typeB=2,r=r.max,r.low=r.min,boot.iter=100)
