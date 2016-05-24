
library(RTDE)


#####
# (1) simulation

N <- 2
n <- 100
x <- RTDE(simu=list(nb=n, marg="ufrechet", cop="indep", replicate=N),
	nbpoint=10:15, alpha=0, omega=1/2)
x	

x$fit$eta
summary(x)
prob(x, 1:4)
summary(prob(x, 1:4))

plot(x, which=1)

y <- RTDE(simu=list(nb=n, marg="ufrechet", cop="indep", replicate=N), 
	nbpoint=10:15, alpha=0:1, omega=1/2)
	
y$fit$eta
	
y
summary(y)
prob(y, 1:4)

plot(y, which=1)


z <- RTDE(simu=list(nb=n, marg="ufrechet", cop="indep", replicate=N), 
	nbpoint=10:15, alpha=0:1, omega=1/2:3, keepdata=TRUE)
z$fit$eta

z
summary(z)
prob(z, 1:4)

plot(z, which=1)


#####
# (2) simulation on multicore

N <- 2
n <- 100
x <- RTDE(simu=list(nb=n, marg="ufrechet", cop="indep", replicate=N),
	nbpoint=10:11, alpha=0, omega=1/2, core=2)
x	
prob(x, 1:4)

x <- RTDE(simu=list(nb=n, marg="ufrechet", cop="indep", replicate=N),
	nbpoint=10:11, alpha=0:1, omega=1/2, core=2)
x	
prob(x, 1:4)
