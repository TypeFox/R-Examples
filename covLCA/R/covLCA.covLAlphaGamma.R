covLCA.covLAlphaGamma <-
function(J,S2,K.j,R,z,y,probs,N,rgivy) #A: Why the difference between the 2 possible computations ?
{
	Cov=matrix(nrow=J*S2*(K.j[1]-1),ncol=J*R*(K.j[1]-1)) #Does not need to be adapted
	
	ind1=0
	for (i1 in 1:J) #A: for each manifest variable m
	{
		for (i2 in 1:S2) #A: for each covariate q
		{
			for (i3 in 1:(K.j[1]-1)) #A: for each category k
			{
				ind1=ind1+1
				ind2=0
				for (i4 in 1:J) #A: for each manifest variable u
				{
					for (i5 in 1:R) #A: for each latent class l
					{
						for (i6 in 1:(K.j[1]-1)) #A: for each category s
						{
							ind2=ind2+1
							Cov[ind1,ind2]= - ((z[,i2]*((y[,i4]==i6)-probs[,i4,i6,i5]))%*%apply((probs[,i1,i3,]*rgivy*(matrix(rep(((1:R)==i5),rep(N,R)),nrow=N,ncol=R)-matrix(rep(rgivy[,i5],R),nrow=N,ncol=R))),1,sum))
							
						}
					}
				}
			}
		}
	}
	return(Cov)
}
