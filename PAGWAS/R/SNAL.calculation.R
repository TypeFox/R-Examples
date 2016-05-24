SNAL.calculation <-function(Y,Phi,s2){
	
	#initialize chain
	z=rep(1,ncol(Phi))

	#first round of iterations
	beta=lars.calculation(x=Phi,y=Y,z=z,s2=s2[1])
	gamma=z^(-1/2)*abs(beta)
	Ggamma=diag(gamma)

	S=Phi%*%Ggamma%*%t(Phi)+diag(length(Y))*s2[1]
	S.inv=matrix.inv.calculation(V=S)

	z=diag(t(Phi)%*%S.inv%*%Phi)
	
	counter=0
		
	repeat{
    		a.check=gamma
    		
    		beta=lars.calculation(x=Phi,y=Y,z=z,s2=s2[1])
    		gamma=z^(-1/2)*abs(beta)
    		Ggamma=diag(gamma)
		S=Phi%*%Ggamma%*%t(Phi)+diag(length(Y))*s2[1]
    		S.inv=matrix.inv.calculation(V=S)
    		z=diag(t(Phi)%*%S.inv%*%Phi)
    		
    		quantity=abs(a.check-gamma)
    		
    		if(length(which(quantity<=1e-6))==length(quantity)){break}
    		counter=counter+1  
    		}

	S.star=Phi%*%diag(gamma)%*%t(Phi)+diag(length(Y))*s2[1]
	S.star.inv=matrix.inv.calculation(V=S.star)
	ARD=as.vector(diag(gamma)%*%t(Phi)%*%S.star.inv%*%Y)
	gamma.star=gamma

	names(gamma.star)=colnames(Phi)
	names(ARD)=colnames(Phi)
	results=list("gamma"=gamma.star,"ARD"=ARD)
	results
}
