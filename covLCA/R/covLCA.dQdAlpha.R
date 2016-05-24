covLCA.dQdAlpha <-
function(rgivy,probs,z,K.j,m,y,S2) #A: checked: same results with 2 different formulae, for both gradient and hessian
{
	N=dim(rgivy)[1]
	R=dim(rgivy)[2]
	Grad=c()
	# Gradbis=c()
	# Gradter=c()
	ind=1
	
	
	for (i1 in 1:S2) #A: for each covariate q
	{
		for (i2 in 1:(K.j[1]-1)) #A: for each category k
		{
			# Formulation 1
			Grad[ind]=(z[,i1]*(y[,m]==i2))%*%apply(rgivy,1,sum) - z[,i1]%*%diag(rgivy%*%t(probs[,m,i2,])) #A: p_imkj
			#Grad[ind]=sum(diag(matrix((y[,m]==i2)*z[,i1],nrow=N,ncol=R)%*%t(rgivy))) - z[,i1]%*%diag(rgivy%*%t(probs[,m,i2,])) #A: p_imkj
			#A: both formulae give the same results

			# Formulation 2 (yields results that are similar, but small differences of order 10^-10
			# Gradbis[ind]=sum(matrix(nrow=N,ncol=R,byrow=FALSE,data=rep(z[,i1],R))*rgivy*(matrix(nrow=N,ncol=R,byrow=FALSE,data=rep((y[,m]==i2),R))-probs[,m,i2,]))
			# Gradter[ind]=(z[,i1]*(y[,m]==i2))%*%apply(rgivy,1,sum) - z[,i1]%*%apply(rgivy*probs[,m,i2,],1,sum)
			
			ind=ind+1
		}
	}
	
	
	Hess=array(dim=c(S2*(K.j[1]-1),S2*(K.j[1]-1))) #Adapted to "S2=1"
	# Hessbis=array(dim=c(S2*(K.j[1]-1),S2*(K.j[1]-1)))
	ind1=0
	ind2=0
	
	for (i1 in 1:S2) #A for each covariate q
	{
		for (i2 in 1:(K.j[1]-1)) #A: for each category k
		{
			ind1=ind1+1
			ind2=0
			for (i3 in 1:S2) #A: for each covariate r
			{
				for (i4 in 1:(K.j[1]-1)) #A: for each category s
				{
					ind2=ind2+1
					Hess[ind1,ind2]=-sum((z[,i1]*z[,i3])%*%(rgivy*probs[,m,i2,]*((i2==i4)-probs[,m,i4,]))) #A: p_imkj
					#Hess[ind1,ind2]=-sum(diag(matrix((z[,i1]*z[,i3]),nrow=N,ncol=R)%*%t(rgivy*probs[,m,i2,]*((i2==i4)-probs[,m,i4,])))) #A: p_imkj
					#A: both formulae give the same results
					
					# Formulation 2 (yields similar results, but differences at 10^-10)
					# Hessbis[ind1,ind2]=-sum(matrix(nrow=N,ncol=R,byrow=FALSE,data=rep(z[,i1],R))*matrix(nrow=N,ncol=R,byrow=FALSE,data=rep(z[,i3],R))*rgivy*probs[,m,i2,]*((i2==i4)-probs[,m,i4,]))
				}
			}
		}
	}
	return(list(grad=Grad,hess=Hess))
}
