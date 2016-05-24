library(boussinesq)

# One-Dimensional Song at al's nonlinear solution 

L <- 1000
x <- seq(from=0,to=L,by=L/100)
t <- c(4,5,20) #  days 

h_sol1 <- beq.song(t=t[1]*3600*24,x=x,s=0.4,h1=1,ks=0.01,nmax=10,alpha=0)
h_sol2 <- beq.song(t=t[2]*3600*24,x=x,s=0.4,h1=1,ks=0.01,nmax=10,alpha=0)
h_sol3 <- beq.song(t=t[3]*3600*24,x=x,s=0.4,h1=1,ks=0.01,nmax=10,alpha=0)
	
	
plot(x,h_sol1,type="l",lty=1,main="Water Surface Elevetion (Song at's solution) ",xlab="x[m]",ylab="h[m]")
lines(x,h_sol2,lty=2)
lines(x,h_sol3,lty=3)
legend("topright",lty=1:3,legend=paste("t=",t,"days",sep=" "))

