covLCA.dQdAlphaGamma <-
function(rgivy,probs,z,K.j,m,S2) #A: rows=gammas, cols=alphas
{
	R=dim(rgivy)[2]
	
	hess=array(dim=c((K.j[1]-1)*R,S2*(K.j[1]-1))) #Adapted to "S2=1"
		hessbis=array(dim=c((K.j[1]-1)*R,S2*(K.j[1]-1)))
	ind1=0
	ind2=0
	
	for (i1 in 1:R) #A: for each LC j
	{
		for (i2 in 1:(K.j[1]-1)) #A: for each category k
		{
			ind1=ind1+1
			ind2=0
			for (i3 in 1:S2)	#A: for each covariate r
			{
				for (i4 in 1:(K.j[1]-1)) #A: for each category s
				{
					ind2=ind2+1
					
					hess[ind1,ind2]=-(z[,i3]*rgivy[,i1])%*%((i2==i4)*probs[,m,i2,i1]-probs[,m,i4,i1]*probs[,m,i2,i1])
					# Formulation 2 (yields very similar results, differences at 10^-12)
					hessbis[ind1,ind2]=-(z[,i3]*rgivy[,i1]*probs[,m,i2,i1])%*%((i2==i4)-probs[,m,i4,i1])
					
				}
			}
		}
	}
	return(hess)
}
