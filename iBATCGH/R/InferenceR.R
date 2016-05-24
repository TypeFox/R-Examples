InferenceR <-
function(listR,G,M,niter,burnin,threshold=0.5){
freTot=RListToVector(listR[(burnin+1):niter],G,M)

FreqMat = NULL
for (i in 1:G){
	FreqMat = rbind(FreqMat,freTot[(1+M*(i-1)):(M*i)])
}
FreqMat=FreqMat/(niter-burnin)

final=matrix(0,G,M)
for (i in 1:G){
for(j in 1:M){
	if(FreqMat[i,j]>threshold){final[i,j]=1}
}
}

out=list(FreqMat=FreqMat,final=final)
return(final)
}
