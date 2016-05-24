# LotkaVolterra.R

# Demo of the simulation functions to analyse the discrete stochastic
# Lotka-Volterra model and compare it to its deterministic
# approximation

# THIS CODE IS VERY SLOW TO RUN!

library(smfsb)
data(spnModels)

# number of samples to plot
samples=20

# number of repeats for computing the mean
repeats=200
#repeats=20

# end point of simulation
endpoint=30

message(paste("Plotting",samples,"individual trajectories"))
stepLV=StepGillespie(LV)
out=simTs(c(x1=50,x2=100),0,endpoint,0.1,stepLV)
plot(out[,2],col="grey",ylim=c(0,800),main="Predator numbers for a stochastic Lotka-Volterra model",ylab="Number of predators")
for (i in 1:samples) {
	out=simTs(c(x1=50,x2=100),0,endpoint,0.1,stepLV)
	lines(out[,2],col="grey")
}

message(paste("Computing the mean of",repeats,"trajectories to overlay"))
cum=simTs(c(x1=50,x2=100),0,endpoint,0.1,stepLV)
names=dimnames(cum)[[2]]
for (i in 2:repeats) {
	out=simTs(c(x1=50,x2=100),0,endpoint,0.1,stepLV)
	cum=cum+out
	message(".",appendLF=FALSE)
}
message("")
cum=cum/repeats
dimnames(cum)[[2]]=names
lines(cum[,2],lwd=2)

message(paste("Computing the mean of",repeats,"trajectories to overlay using C function"))
cum=simTs(c(x1=50,x2=100),0,endpoint,0.1,stepLVc)
names=dimnames(cum)[[2]]
for (i in 2:repeats) {
	out=simTs(c(x1=50,x2=100),0,endpoint,0.1,stepLVc)
	cum=cum+out
	message(".",appendLF=FALSE)
}
message("")
cum=cum/repeats
dimnames(cum)[[2]]=names
lines(cum[,2],lwd=2,col=2)

message(paste("Computing the mean of",repeats,"CLE trajectories to overlay as a dashed line"))
stepLVCLE=StepGillespie(LV)
cum=simTs(c(x1=50,x2=100),0,endpoint,0.1,stepLVCLE)
names=dimnames(cum)[[2]]
for (i in 2:repeats) {
	out=simTs(c(x1=50,x2=100),0,endpoint,0.1,stepLVCLE)
	cum=cum+out
	message(".",appendLF=FALSE)
}
message("")
cum=cum/repeats
dimnames(cum)[[2]]=names
lines(cum[,2],lty=2,lwd=2)

message("Overlaying the deterministic solution as a dotted line")
stepLVEuler=StepEulerSPN(LV,dt=0.0001)
out=simTs(c(x1=50,x2=100),0,endpoint,0.1,stepLVEuler)
lines(out[,2],lty=3,lwd=2)



# eof

