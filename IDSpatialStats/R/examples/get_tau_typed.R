data(DengueSimulationR02)

r.max<-seq(20,1000,20)
r.min<-seq(0,980,20)
r.mid<-(r.max+r.min)/2

#Lets see if there's a difference in spatial dependence by time case occurs
type<-2-(DengueSimR02[,"time"]<120)
tmp<-cbind(DengueSimR02,type=type)

typed.tau<-get.tau.typed(tmp,typeA=1,typeB=2,r=r.max,r.low=r.min,comparison.type = "independent")

plot(r.mid,typed.tau,log="y",cex.axis=1.25,
     xlab="Distance (m)",ylab="Tau",cex.main=0.9,lwd=2,type="l")
abline(h=1,lty=2)
