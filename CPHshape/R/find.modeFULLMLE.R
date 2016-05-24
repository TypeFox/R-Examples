find.modeFULLMLE <-
function(x, z, delta, mode, beta0=rep(1, length(as.matrix(z)[1,])), max.loop=250, eps=1e-5, eps.beta=1e-15, print=FALSE){

	u			<-	order(x)
	x			<-	x[u]
	delta		<-	delta[u]
	z			<-	z[u,]
                                
	beta.new	<-	beta0
	mle			<-	find.modeMLE(x, w=as.vector(exp(z%*%as.matrix(beta.new))), mode, delta)
	phi.new		<-	mle$phi-sum(delta*as.vector((z%*%as.matrix(beta.new))))
	H.new		<-	mle$H     

	if(print==TRUE) { cat("iter=i", "phi[i]", "|phi[i]-phi[i-1]|", "beta(s)", "\n") }
	if(print==TRUE){cat(0, phi.new, "NA", beta.new, "\n")}
                                                                        
	for (i in 1:max.loop){
                
		phi.old		<-	phi.new
		beta.old	<-	beta.new
		H.old		<-	H.new
                
		beta.new	<-	find.beta(x, delta, z, H.old, max.loop=max.loop, eps=eps.beta)
		mle			<-	find.modeMLE(x, w=as.vector(exp(z%*%as.matrix(beta.new))), mode, delta)
		H.new		<-	mle$H
		phi.new		<-	mle$phi-sum(delta*as.vector((z%*%as.matrix(beta.new))))
                
		if(print==TRUE){cat(i, phi.new, abs(phi.new-phi.old), beta.new, "\n")}
		if(abs(phi.new-phi.old) < eps) {break}

		}  # end loop
                
	return(list(beta=beta.new, h.range=mle$ranges, h.val=mle$mle, phi=phi.new, H=mle$H, mode=mode))
	}  # find.modeFULLMLE

