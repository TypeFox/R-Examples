covLCA.identifiability <-
function(J,K.j,R,x,z,npar,alph,bet,gamm)
{
	#A: condition 1: the number of unique model parameters cannot exceed the number of independent pieces of observed information.
	if (npar > (prod(K.j)-1))
	{
		cat("\n ALERT: the model is not locally identifiable (the number of parameters exceeds the number of independent pieces of observed information) \n \n")
	}
	
	#A: condition 2: free model parameters (alphas, gammas, betas) and covariate values (x,z) are all finite
	if ( any(abs(alph)==Inf) | any(abs(gamm)==Inf) | any(abs(bet)==Inf))
	{
		cat("\n ALERT: the model is not locally identifiable (at least one parameter is not finite) \n \n")
	}
	
	if ( any(abs(x)==Inf) | any(abs(z)==Inf) )
	{
		cat("\n ALERT: the model is not locally identifiable (at least one covariate value is not finite) \n \n")
	}
	
	#A: condition 3: vectors tau are linearly independent
	library(Matrix)
	gam=array(gamm,dim=c(J,K.j[1]-1,R)) #A: gamma_mkj
	
	tauMat=matrix(nrow=prod(K.j)-1,ncol=R) #A: rows=response patterns, cols=LC
	gamma2=array(dim=c(J,K.j[1],R))

	datacell=matrix(nrow=prod(K.j),ncol=J)
	for (m in 1:J) #A: for each manifest variable
	{
		datacell[,(J-m+1)]=rep(seq(1,K.j[1]),rep(K.j[1]^(m-1),K.j[1]))
	}
	datac=datacell[1:(dim(datacell)[1]-1),] #A: don't consider the last (reference) profile
	
	for (j in 1:R) #A: for each latent class
	{
		#gamma2[,,j]=cbind(gamm[,,j],0) #A: gamma_mKj=0 forall m,j
		gamma2[,,j]=cbind(gam[,,j],0) #A: gamma_mKj=0 forall m,j
		tauMat[,j]=rep(1,prod(K.j)-1)
		for (m in 1:J) #A: for each manifest variable
		{
			tauMat[,j]=tauMat[,j]*( exp(gamma2[m,datac[,m],j])/(1+sum(exp(gam[m,,j]))) )#A: num: gamma2, denom: gam
		}
	}
	tauMat.eigen=eigen(t(tauMat)%*%tauMat,only.values=TRUE)$values
	tau.invcond=rcond(tauMat)
	
	#if(any(tauMat.eigen < 0.0000001)) 
	#{
	#	cat("\n ALERT: the model is not locally identifiable (tau vectors are not linearly independent) \n \n")
	#}
	
	if (tau.invcond< 0.00000001)
	{
		cat("\n ALERT: the model is not locally identifiable (tau vectors are not linearly independent \n \n")
	}
	
	#A: condition 4: both design matrices x and z have full column rank
	xmat=as.matrix(x)	
	x.eigen=eigen(t(xmat)%*%xmat,only.values=TRUE)$values
	x.invcond=rcond(xmat)
	z.eigen=eigen(t(cbind(1,as.matrix(z)))%*%cbind(1,as.matrix(z)),only.values=TRUE)$values
	z.invcond=rcond(cbind(1,as.matrix(z))) #added cbind(1, ) on 30JUL
	#if(any(x.eigen < 0.0000001)) 
	#{
	#	cat("\n ALERT: the model is not locally identifiable (matrix x has not full column rank) \n \n")
	#}
	
	#if(any(z.eigen < 0.0000001)) 
	#{
	#	cat("\n ALERT: the model is not locally identifiable (matrix z has not full column rank) \n \n")
	#}
		
	if(x.invcond < 0.00000001)
	{
		cat("\n ALERT: the model is not locally identifiable (matrix x has not full column rank) \n \n")
	}
	if(z.invcond < 0.00000001) 
	{
		cat("\n ALERT: the model is not locally identifiable (matrix z has not full column rank) \n \n")
	}	
	
	
	return(list(tau.eigen=tauMat.eigen,x.eigen=x.eigen,z.eigen=z.eigen,tau.invCond=tau.invcond,x.invCond=x.invcond,z.invCond=z.invcond))	
}
