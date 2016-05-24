covLCA.covLBetaBeta <-
function(rgivy,x,S1,R)
{
	Cov=matrix(nrow=S1*(R-1),ncol=S1*(R-1))
	ind1=0
	for (i1 in 1:(R-1)) #A: for each latent class j
	{
		for (i2 in 1:S1) #A: for each covariate p
		{
			ind1=ind1+1
			ind2=0
			for (i3 in 1:(R-1)) #A: for each latent class l
			{
				for (i4 in 1:S1) #A: for each covariate u
				{
					ind2=ind2+1
					Cov[ind1,ind2]=(x[,i2]*x[,i4]*rgivy[,i1])%*%((i1==i3)-rgivy[,i3])					
				}
			}
		}
	}
	return(Cov)
}
