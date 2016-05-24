covLCA.covLBetaAlpha <-
function(S1,R,J,S2,K.j,x,z,probs,rgivy,N)
{
	Cov=matrix(nrow=S1*(R-1),ncol=J*S2*(K.j[1]-1)) #Does not need to be adapted (bc J>1)
	
	ind1=0
	for (i1 in 1:(R-1)) #A: for each LC j
	{
		for (i2 in 1:S1) #A: for each covariate p
		{
			ind1=ind1+1
			ind2=0
			for (i3 in 1:J) #A: for each manifest variable u
			{
				for (i4 in 1:S2) #A: for each covariate r
				{
					for (i5 in 1:(K.j[1]-1)) #A: for each category s
					{
						ind2=ind2+1
						Cov[ind1,ind2]= - ((x[,i2]*z[,i4]) %*% apply((probs[,i3,i5,]*matrix(rep(rgivy[,i1],R),nrow=N,ncol=R)*(matrix(rep((i1==(1:R)),rep(N,R)),nrow=N,ncol=R)-rgivy)),1,sum))
						
						
					}
				}
			}
		}
	}
	return(Cov)
}
