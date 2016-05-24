# ZeBook # Francois Brun
# version 2013-05-23
################################ MAIN program ##################################
library(ZeBook)
################################################################################
# main program
a=0.1
Y0=1.0
duration=40
dt=2

# a) comparison between analytical solution
# comparison between analytical solution and the numeric Euler integration
x=seq(0,duration,by=0.1)
y=Y0*exp(a*x)
plot(x,y,type="l",lwd=3, xlab="time (day)",ylab="Y (density of population)")
lines(exponential.model(a,Y0,duration,dt),lty=2,lwd=3,col="black")
legend("topleft",legend=c("Analytical Solution","Numerical Solution"),lwd=c(2,2),lty=c(1,2),col=c("black","black"),cex=0.75)

dev.new()
# comparison between analytical solution and the numeric improved Euler integration
x=seq(0,duration,by=0.1)
y=Y0*exp(a*x)
plot(x,y,type="l",lwd=3, xlab="time (day)",ylab="Y (density of population)")
lines(exponential.model.ie(a,Y0,duration,dt),lty=2,lwd=3,col="black")
legend("topleft",legend=c("Analytical Solution","Numerical Solution"),lwd=c(2,2),lty=c(1,2),col=c("black","black"),cex=0.75)

# b) influence of dt on Euler integration : error of integration
list_dt=c(2,1,0.5,0.1,0.01)
sim_dt = data.frame(time=seq(0,duration,by=2))
for(dt in list_dt){
    sim = exponential.model(a,Y0,duration=duration,dt=dt)
    sim_dt = cbind(sim_dt, Y=sim[sim$time%in%sim_dt$time ,"Y"])
 }
names(sim_dt)[-1]=paste("Y_dt", list_dt,sep="")

dev.new()
plot(sim_dt$time, sim_dt$Y_dt0.01,type="l",lwd=2,xlim=c(0,max(sim_dt$time)+10), ylim=c(0,max(sim_dt$Y_dt0.01)+0.2),xlab="time (day)",ylab="Y (density of population)")
lines(sim_dt$time, sim_dt$Y_dt0.1, lty=2,col="darkgrey")
lines(sim_dt$time, sim_dt$Y_dt0.5, lty=2,col="darkgrey")
lines(sim_dt$time, sim_dt$Y_dt1, lty=2,col="darkgrey")
lines(sim_dt$time, sim_dt$Y_dt2, lty=2,col="darkgrey")
text(max(sim_dt$time),max(sim_dt$Y_dt0.01)+0.1,"dt=0.01",pos=4,cex=0.8 )
text(max(sim_dt$time),max(sim_dt$Y_dt0.1),"dt=0.1",pos=4,cex=0.8,col="darkgrey")
text(max(sim_dt$time),max(sim_dt$Y_dt0.5),"dt=0.5",pos=4,cex=0.8,col="darkgrey")
text(max(sim_dt$time),max(sim_dt$Y_dt1),"dt=1",pos=4,cex=0.8,col="darkgrey")
text(max(sim_dt$time),max(sim_dt$Y_dt2),"dt=2",pos=4,cex=0.8,col="darkgrey")

# c) influence of dt on Euler integration : plot of error of integration
list_dt=c(0.01,0.1,0.5,1,2,2.5,4.0)
sim_dt = data.frame(time=seq(0,duration,by=0.5))
for(dt in list_dt){
    sim = exponential.model(a,Y0,duration=duration,dt=dt)
    sim_dt = merge(sim_dt, sim, by="time",all.x=TRUE)
 }
names(sim_dt)[-1]=paste("Y_dt", list_dt,sep="")
dev.new()

# compute Yexact - YEuler at time=40
sim_dt_t40= sim_dt[sim_dt$time==40,][-1]
error= Y0*exp(a*40) - sim_dt_t40
plot(list_dt, error,type="l",lwd=2,xlim=c(0,4), ylim=range(error,na.rm=TRUE),xlab="dt",ylab="Error of integration")
# end of file
