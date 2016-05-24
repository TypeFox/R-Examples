covLCA.covLAlphaAlpha <-
function(J,S2,K.j,N,R,probs,rgivy,z)
{
	Cov=matrix(nrow=J*S2*(K.j[1]-1),ncol=J*S2*(K.j[1]-1)) #Does not need to be adapted
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
					for (i5 in 1:S2) #A: for each covariate r
					{
						for (i6 in 1:(K.j[1]-1)) #A: for each category s
						{
							ind2=ind2+1
							res=rep(0,N)
							for (j in 1:R) #A: for each latent class j
							{
								for (l in 1:R) #A: for each latent class l
								{
									res=res+probs[,i1,i3,j]*probs[,i4,i6,l]*rgivy[,j]*((j==l)-rgivy[,l])
								}
							}
							Cov[ind1,ind2]=(z[,i2]*z[,i5])%*%res
							
						}
					}
				}
			}
		}
	}
	return(Cov)
}
