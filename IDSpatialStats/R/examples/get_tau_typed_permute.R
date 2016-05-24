data(DengueSimulationR02)

r.max<-seq(20,1000,20)
r.min<-seq(0,980,20)
r.mid<-(r.max+r.min)/2

#Lets see if there's a difference in spatial dependence by time case occurs
type<-2-(DengueSimR02[,"time"]<120)
tmp<-cbind(DengueSimR02,type=type)

typed.tau<-get.tau.typed(tmp,typeA=1,typeB=2,r=r.max,r.low=r.min,comparison.type = "independent")
typed.tau.type.null<-get.tau.typed.permute(tmp,typeA=1,typeB=2,r=r.max,r.low=r.min,permutations=100,comparison.type = "independent")

null.ci<-apply(typed.tau.type.null,2,quantile,probs=c(0.025,0.975))

plot(r.mid,typed.tau,ylim=c(0.3,4),log="y",cex.axis=1.25,,xlab="Distance (m)",ylab="Tau",cex.main=0.9,lwd=2,type="n")
abline(h=1,lty=1)
lines(r.mid,typed.tau,pch=20,col=1,lwd=3)
lines(r.mid, null.ci[1,] , lty=2)
lines(r.mid, null.ci[2,] , lty=2)

