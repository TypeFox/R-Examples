covLCA.dQdBeta <-
function(rgivy,prior,x) #A: rgivy: posterior probs, a matrix where rows=indiv, cols=LC ; prior: rows=indiv, cols=LC
{	
	
	Grad=c()
	ind=1
	for (i2 in 1:(dim(prior)[2]-1)) #A: for each LC (except the last one)
	{
		for (i1 in 1:dim(x)[2]) #A: for each covariate
		{
			Grad[ind]=x[,i1]%*%(rgivy[,i2]-prior[,i2])
			ind=ind+1
		}
	} #A: grad contains derivatives wrt beta_jp where jp= 11,12,13,...,1P, 21,22,...,2P;... J1,...,JP
	
	Hess=matrix(nrow=dim(x)[2]*(dim(prior)[2]-1),ncol=dim(x)[2]*(dim(prior)[2]-1))
	
	ind1=0
	ind2=0
	
	for (i2 in 1:(dim(prior)[2]-1)) #A: for each LC j(rows) (except the last one)
	{
		for (i1 in 1:dim(x)[2]) #A: for each covariate p (rows)
		{
			ind1=ind1+1
			ind2=0
			for (i4 in 1:(dim(prior)[2]-1)) #A: for each LC l (cols)
			{
				for (i3 in 1:dim(x)[2]) #A: for each covariate u(cols)
				{
					ind2=ind2+1
					Hess[ind1,ind2]=-(x[,i1]*x[,i3])%*% ( prior[,i2]*((i2==i4)-prior[,i4]) )
				}
			}
		}
	}
	return(list(grad=Grad,hess=Hess))
}
