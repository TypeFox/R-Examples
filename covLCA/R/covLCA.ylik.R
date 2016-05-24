covLCA.ylik <-
function(probs,y,K.j) 
{
	J=dim(y)[2]
	N=dim(y)[1]
	classes=dim(probs)[4]
	y.imk=array(0,dim=c(N,J,K.j[1])) #A: a 3D array wich will  contain y_imk
	
	
	for (i1 in 1:J) #A: for each manifest variable
	{
		for (i2 in 1:K.j[1]) #A: for each category of the manifest variable
		{
			y.imk[,i1,i2]=(y[,i1]==i2)
		}
	}
	y.imk.rep=array(y.imk,dim=c(N,J,K.j[1],classes)) #A:this 4D-array contains, for each LC (4th dimension), the 3D-array y.imk
	mat=probs^y.imk.rep #A: mat is a 4D-array which contains p_imkj^(y_imk)
	
	lik=matrix(1,nrow=dim(y)[1],ncol=classes) #A: rows=indiv,cols=LC
	for (i3 in 1:classes) #A: for each LC
	{
		lik[,i3]=apply(mat[,,,i3],1,prod) #A: make product for each individual, over manifest variables and LC
	}
	return(lik)
	
}
