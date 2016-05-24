fun.ntitr0 <-
function(u0,m,b,oz,od,Z,n){

	## ntitr0.m
	#------------------------------------#

	#data10 <- fun.oldp2(b,m,Data)
	#s <- data10$s
	#ru <- data10$ru 
	#u <- data10$u

	ru <- u0 %*% matrix(1,nrow = 1, ncol = m)
	gama1 <- exp(-oz %*% b[1,])
    gama2 <- exp(-oz %*% b[2,])

	a <- gama1+gama2*ru

	s <- t(od)%*%log(a) + colSums(log(a/gama1)/gama2)

	ow <- matrix(1,nrow = n, ncol = 1)

	u1 <- t(oz*od)%*%(-gama1/a)+t(Z)%*%(ru/a)
    u2 <- -t(oz*od)%*%(gama2*ru/a)-t(Z)%*%(ru/a)+t(oz)%*%(log(a/gama1)/gama2)

	u <- rbind(u1,u2)
	#------------------------------------#

#-----------------------------------------------------------------#
## Output Resuts 
#-----------------------------------------------------------------#
    output<- list(u=u,s=s,ru=ru)
    return(output)
#-----------------------------------------------------------------#
}
