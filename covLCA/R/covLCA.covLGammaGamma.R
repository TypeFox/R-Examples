covLCA.covLGammaGamma <-
function(J,R,K.j,probs,rgivy,y)
{
	Cov=matrix(nrow=J*R*(K.j[1]-1),ncol=J*R*(K.j[1]-1))
	ind1=0
	for (i1 in 1:J) #A: for each manifest variable m
	{
		for (i2 in 1:R) #A: for each latent class j
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
							Cov[ind1,ind2]=(((y[,i1]==i3)-probs[,i1,i3,i2])*((y[,i4]==i6)-probs[,i4,i6,i5]))%*%(rgivy[,i2]*((i2==i5)-rgivy[,i5]))
						}
					}
				}
			}
		}
	}
	return(Cov)
}
