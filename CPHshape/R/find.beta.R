find.beta <-
function(x, delta, z, H, beta0=rep(1,length(as.matrix(z)[1,])), max.loop=250, eps=10^-10, print=FALSE){

	z		<-	as.matrix(z)
	beta	<-	as.matrix(beta0)
	phi.new	<-	sum(H*exp(z%*%as.matrix(beta))-delta*(z%*%as.matrix(beta)))

	for (i in 1:max.loop){

		grad	<-	t(z)%*%(H*exp(z%*%beta)-delta)
		hess	<-	t(z)%*%diag(as.vector(H*exp(z%*%beta)))%*%z
		beta	<-	beta-solve(hess)%*%grad
		phi.old <-	phi.new
		phi.new	<-	sum(H*exp(z%*%beta)-delta*(z%*%beta))
                                
		if(print==TRUE) { cat(i, phi.old, abs(phi.new-phi.old), grad, as.vector(beta), "\n") }
		if(sqrt(sum(grad^2))< eps){break}
		}

	return(as.vector(beta))
	}

