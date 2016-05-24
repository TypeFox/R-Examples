covLCA.updateCond <-
function(a,g,z,R,J,K.j,S2,N) 
{
	gama <- array(g,dim=c(J,K.j[1]-1,R)) #A: gamma_mkj
	alpha <- array(a,dim=c(J,K.j[1]-1,S2)) #A: alpha_mkp
	
	p=array(dim=c(N,J,K.j[1],R)) #A: p_imkj
	# pbis=array(dim=c(N,J,K.j[1],R))
		
	for (i2 in 1:J) #A: for each manifest variable
	{
		for (i3 in 1:R) #A: For each LC 
		{

			# 1st formulation
			#p[,i2,,i3]=exp(matrix(rep(c(gama[i2,,i3],0),rep(N,K.j[1])),nrow=N,ncol=K.j[1])+z%*%t(rbind(alpha[i2,,],0)))/apply(exp(matrix(rep(c(gama[i2,,i3],0),rep(N,K.j[1])),nrow=N,ncol=K.j[1])+z%*%t(rbind(alpha[i2,,],0))),1,sum)
			if(dim(alpha)[3]!=1) {p[,i2,,i3]=exp(matrix(rep(c(gama[i2,,i3],0),rep(N,K.j[1])),nrow=N,ncol=K.j[1])+z%*%t(rbind(alpha[i2,,],0)))/apply(exp(matrix(rep(c(gama[i2,,i3],0),rep(N,K.j[1])),nrow=N,ncol=K.j[1])+z%*%t(rbind(alpha[i2,,],0))),1,sum)}
			if(dim(alpha)[3]==1) {p[,i2,,i3]=exp(matrix(rep(c(gama[i2,,i3],0),rep(N,K.j[1])),nrow=N,ncol=K.j[1])+z%*%array(c(alpha[i2,,],0),dim=K.j[1]))/apply(exp(matrix(rep(c(gama[i2,,i3],0),rep(N,K.j[1])),nrow=N,ncol=K.j[1])+z%*%array(c(alpha[i2,,],0),dim=K.j[1])),1,sum)}

			
			#A: repeat vector gamma_m.j, of length K, because to each column of the following matrix, we must add the corresponding element of this vector
			#p[i1,i2,i3,K.j]=1/(1+sum(exp(gamm[i2,,i3]+z[i1,]%*%t(alph[i2,,])))) 
			#A: last category: gamma and alphas = 0
			
			# 2nd formulation (yields the same results as formulation 1)
			# temp=exp(cbind(matrix(nrow=N,ncol=K.j[1]-1,byrow=TRUE,data=rep(gama[i2,,i3],N)),0)+z%*%cbind(t(alpha[i2,,]),0))
			# temp2=apply(temp,1,sum)
			# pbis[,i2,,i3]=temp/temp2
		}
	}
	
	return(p)	
}
