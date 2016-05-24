# ExtrinsicNoise.R
# Immigration-death diffusion approx with death rate a CIR process

require(smfsb)

theta=c(lambda=2,alpha=1,mu=0.1,sigma=0.2)

myDrift <- function(x,t,th=theta)
     {
             with(as.list(c(x,th)),{
                     c( lambda - x*y ,
                           alpha*(mu-y) )
             })
     }

myDiffusion <- function(x,t,th=theta)
     {
             with(as.list(c(x,th)),{
                     matrix(c( sqrt(lambda + x*y) , 0,
                           0, sigma*sqrt(y) ),ncol=2,nrow=2,byrow=TRUE)
             })
     }

diffusionE <- function(x,t,th=theta)
     {
             with(as.list(c(x,th)),{
                     matrix(c( 0 , 0,
                           0, sigma*sqrt(y) ),ncol=2,nrow=2,byrow=TRUE)
             })
     }

diffusionI <- function(x,t,th=theta)
     {
             with(as.list(c(x,th)),{
                     matrix(c( sqrt(lambda + x*y) , 0,
                           0, 0 ),ncol=2,nrow=2,byrow=TRUE)
             })
     }

dt=0.001
op=par(mfrow=c(2,2))
stepProc=StepEuler(myDrift,dt=dt)
stepProcE=StepSDE(myDrift,diffusionE,dt=dt)
stepProcI=StepSDE(myDrift,diffusionI,dt=dt)
stepProcIE=StepSDE(myDrift,myDiffusion,dt=dt)
for (stepFun in c(stepProc,stepProcE,stepProcI,stepProcIE)) {
	out=simTs(c(x=1,y=0.1),0,30,0.01,stepFun)
	plot(out[,1],ylim=c(0,30),ylab="x(t)",lwd=2)
	}
par(op)

# eof

