covLCA.covLBetaGamma <-
function(S1,R,J,K.j,x,y,probs,rgivy)
{
	Cov=matrix(nrow=S1*(R-1),ncol=J*R*(K.j[1]-1))
	ind1=0
	for (i1 in 1:(R-1)) #A: for each LC j
	{
		for (i2 in 1:S1) #A: for each covariate p
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
						Cov[ind1,ind2]= (x[,i2]*((y[,i4]==i6)-probs[,i4,i6,i5]))%*%(rgivy[,i1]*((i1==i5)-rgivy[,i5])) #A: error corrected
					}
				}
			}
		}
	}
	return(Cov)
}
