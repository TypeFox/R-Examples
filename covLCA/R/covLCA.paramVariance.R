covLCA.paramVariance <-
function(prior,probs,rgivy,R,S1,S2,J,K.j,x,y,z,N)
{
	
	ret=ginv(-covLCA.varianceFirst(prior,probs,rgivy,R,S1,J,K.j,S2,x,y,z)-covLCA.varianceSecond(rgivy,x,S1,R,J,K.j,probs,y,S2,N,z))
	return(ret)
}
