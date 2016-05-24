covLCA.postClass <-
function(prior,probs,y,K.j) #A: prior : rows=indiv, cols=LC ; ylik(vp,yp): a matrix where rows=indiv, cols=LC
{
	N=dim(y)[1]
	J=dim(y)[2]
	R=dim(probs)[4]
	postProbs=matrix(nrow=dim(y)[1],ncol=R)
	# postProbsbis=matrix(nrow=dim(y)[1],ncol=R)

	# Formulation 1
	numer=prior*covLCA.ylik(probs,y,K.j)
	denom=apply(numer,1,sum) #A: for each individual, sum over the LC
	postProbs=numer/matrix(denom,ncol=R,nrow=N) #A: matrix where rows=individuals, columns=latent classes. The denom is the same for each indiv(row)
	
	# Formulation 2 (yields same results as 1st formulation
	# temp3=covLCA.ylik(probs,y,K.j)*prior
	# temp4=apply(temp3,1,sum)
	# postProbsbis=temp3/temp4
	
	
	return(postProbs)
}
