library(boussinesq)

# One dimensinal linearized Boussinesq Equation (ground-water)

L <- 1000
x <- seq(from=0,to=L,by=L/100)
t <- 4 # 4 days 
h_sol0 <- beq.lin(x=x,t=t*24*3600,h1=2,h2=1,ks=0.01,L=L,s=0.4,big=100,p=0.0)
h_solp <- beq.lin(x=x,t=t*24*3600,h1=2,h2=1,ks=0.01,L=L,s=0.4,big=100,p=0.5)
h_sol1 <- beq.lin(x=x,t=t*24*3600,h1=2,h2=1,ks=0.01,L=L,s=0.4,big=100,p=1.0)

plot(x,h_sol0,type="l",lty=1,main=paste("Water Surface Elevetion after",t,"days",sep=" "),xlab="x[m]",ylab="h[m]")
lines(x,h_solp,lty=2)
lines(x,h_sol1,lty=3)
legend("topright",lty=1:3,legend=c("p=0","p=0.5","p=1"))

