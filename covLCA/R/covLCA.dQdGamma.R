covLCA.dQdGamma <-
function(rgivy,probs,y,K.j,m)
{
	R=dim(rgivy)[2] #A: nb of LC
	
	Grad=c()
	ind=1

	for (i2 in 1:R) #A: for each LC 
	{
		for (i3 in 1:(K.j[1]-1)) #A: for each category (except the last one)
		{
			#Grad[ind]= rgivy[,i2]%*%((y[,m]==i3)-probs[,m,i3,i2]) #A: p_imkj
			Grad[ind]= sum( rgivy[,i2]*((y[,m]==i3)-probs[,m,i3,i2] ))
			#A: no difference
			ind=ind+1
		}
	}
		
	Hess=matrix(nrow=(K.j[1]-1)*R,ncol=(K.j[1]-1)*R)
	ind1=0
	ind2=0
	
	for (i1 in 1:R) #A: for each LC j
	{
		for (i2 in 1:(K.j[1]-1)) #A: for each category k
		{
			ind1=ind1+1
			ind2=0
			for (i3 in 1:R) #A: for each LC l
			{
				for (i4 in 1:(K.j[1]-1)) #A: for each category s
				{
					ind2=ind2+1
					Hess[ind1,ind2]= - (rgivy[,i1] * (i1==i3) * probs[,m,i2,i1])%*% ((i2==i4)-probs[,m,i4,i3])
				}
			}
		}
	}
	return(list(grad=Grad,hess=Hess))
}
