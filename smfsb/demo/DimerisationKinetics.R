# DimerisationKinetics.R
# Demo of the dimerisation kinetics model from chapter 7

require(smfsb)
data(spnModels)

stepDimer=StepGillespie(Dimer)
op=par(mfrow=c(2,2))
out=simTs(c(x1=301,x2=0),0,10,0.01,stepDimer)
plot(out,plot.type="single")
plot(out[,1],ylim=c(0,300))
for (i in 1:19) {
	out=simTs(c(x1=301,x2=0),0,10,0.01,stepDimer)
	lines(out[,1])
}
hist(simSample(1000,c(x1=301,x2=0),0,10,stepDimer)[,1],30)
par(op)


# eof

