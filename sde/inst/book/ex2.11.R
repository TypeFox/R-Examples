# ex2.11.R
S0 <- 1
theta <- c(1, 0.5)
euler <- NULL
milstein <- NULL
for(i in 1:14){
	n <- length(W[[i]])
	dt <- 1/n
	sdt <- sqrt(dt)
	E <- numeric(n)
	E[1] <- S0 
	M <- numeric(n)
	M[1] <- S0 
	for(j in 2:n){
		Z <- W[[i]][j]-W[[i]][j-1]
		E[j] <-  E[j-1] * (1 + theta[1] * dt + theta[2] * Z)  
		M[j] <-  M[j-1] * (1 + (theta[1] - 0.5* theta[2]^2) * dt + 
		 theta[2] * Z + 0.5 * theta[2]^2 * Z^2)  
	}
	cat(paste(E[n],M[n],"\n"))
	euler <- c(euler, E[n])
	milstein <- c(milstein, M[n])
}
plot(1:14,euler,type="l",main="Milstein vs Euler",
 xlab=expression(log[2](N)), ylab="S(T)")
lines(1:14,milstein,lty=2)
abline(h=exp(theta[1]-0.5*theta[2]^2),lty=3)
