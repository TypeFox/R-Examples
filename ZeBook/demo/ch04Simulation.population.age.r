# ZeBook # Francois Brun
# version 2013-03-07
################################ MAIN program ##################################
library(ZeBook)
# 1) simulation
system.time(sim <- population.age.model(rb=3.5,mE=0.017,rE=0.172,m1=0.060,r12=0.217,m2=0.032,r23=0.313,m3=0.022,r34=0.222,m4=0.020,r4P=0.135,mP=0.020,rPA=0.099,mA=0.027,iA=0,duration=40,dt=0.01))

# 2) values simulated at different time:
sim[sim$time==0,]
sim[sim$time==10,]
sim[sim$time==40,]

# 3) Graphical representation
graph.param = data.frame("V"=c("E","L1","L2","L3","L4","P","A"),
"lty"=c(2,2,2,2,3,4,1), "lwd"=c(2,1,1,1,1,2,3))
plot(c(0,max(sim$time)), c(0,max(sim[,-1])),type="n",xlab="time (day)",ylab="population density")
null=sapply(c("E","L1","L2","L3","L4","P","A"),function(v) lines(sim$time,sim[,v],
lty=graph.param[graph.param$V==v,"lty"],lwd=graph.param[graph.param$V==v,"lwd"]))
legend("topright", legend=graph.param$V, lty = graph.param$lty, lwd = graph.param$lwd, cex=0.75)


# 4) influence of dt on Euler integration : error of integration
list_dt=c(2,1,0.5,0.1,0.01)
sim_dt = data.frame(time=seq(0,40,by=2))
for(dt in list_dt){
    sim = population.age.model(rb=3.5,mE=0.017,rE=0.172,m1=0.060,r12=0.217,m2=0.032,r23=0.313,m3=0.022,r34=0.222,m4=0.020,r4P=0.135,mP=0.020,rPA=0.099,mA=0.027,duration=max(sim_dt$time),dt=dt)
    sim_dt = cbind(sim_dt, A=sim[sim$time%in%sim_dt$time ,"A"])
 }
names(sim_dt)[-1]=paste("A_dt", list_dt,sep="")

dev.new()
plot(sim_dt$time, sim_dt$A_dt0.01,type="l",lwd=2,xlim=c(0,max(sim_dt$time)+10), ylim=c(0,max(sim_dt$A_dt0.01)+0.2),xlab="time (day)",ylab="A (density of population)")
lines(sim_dt$time, sim_dt$A_dt0.1, lty=2,col="grey")
lines(sim_dt$time, sim_dt$A_dt0.5, lty=2,col="grey")
lines(sim_dt$time, sim_dt$A_dt1, lty=2,col="grey")
lines(sim_dt$time, sim_dt$A_dt2, lty=2,col="grey")
text(max(sim_dt$time),max(sim_dt$A_dt0.01)+0.1,"dt=0.01",pos=4,cex=0.8 )
text(max(sim_dt$time),max(sim_dt$A_dt0.1),"dt=0.1",pos=4,cex=0.8,col="grey")
text(max(sim_dt$time),max(sim_dt$A_dt0.5),"dt=0.5",pos=4,cex=0.8,col="grey")
text(max(sim_dt$time),max(sim_dt$A_dt1),"dt=1",pos=4,cex=0.8,col="grey")
text(max(sim_dt$time),max(sim_dt$A_dt2),"dt=2",pos=4,cex=0.8,col="grey")

# 5) comparison with the matrix writen model : performance
system.time(simM <- population.age.matrix.model(rb=3.5,mE=0.017,rE=0.172,m1=0.060,r12=0.217,m2=0.032,r23=0.313,m3=0.022,r34=0.222,m4=0.020,r4P=0.135,mP=0.020,rPA=0.099,mA=0.027,iA=0,duration=150,dt=0.1)) 

system.time(sim <- population.age.model(rb=3.5,mE=0.017,rE=0.172,m1=0.060,r12=0.217,m2=0.032,r23=0.313,m3=0.022,r34=0.222,m4=0.020,r4P=0.135,mP=0.020,rPA=0.099,mA=0.027,iA=0,duration=150,dt=0.1))

all(sim==simM)

# 6) simulation using the model written with ode 
library(deSolve)

system.time(simode <- population.age.model.ode(rb=3.5,mE=0.017,rE=0.172,m1=0.060,r12=0.217,m2=0.032,r23=0.313,m3=0.022,r34=0.222,m4=0.020,r4P=0.135,mP=0.020,rPA=0.099,mA=0.027,iA=0,duration=40,dt=1, method="rk4"))
# end of file
