InferenceXi <-
function(listXi,niter,burnin){
FreqMat1=listXi[[1]]
FreqMat3=listXi[[2]]
FreqMat4=listXi[[3]]

s=dim(FreqMat1)[[1]]
m=dim(FreqMat1)[[2]]

FreqMat2=matrix(niter-burnin,nrow=s,ncol=m)
FreqMat2=FreqMat2-FreqMat1-FreqMat3-FreqMat4

FreqMat=matrix(0,nrow=s,ncol=m)
for(i in 1:s){
  for(j in 1:m){
    if(max(c(FreqMat1[i,j],FreqMat2[i,j],FreqMat3[i,j],FreqMat4[i,j]))==FreqMat1[i,j]){FreqMat[i,j]=1}
    if(max(c(FreqMat1[i,j],FreqMat2[i,j],FreqMat3[i,j],FreqMat4[i,j]))==FreqMat2[i,j]){FreqMat[i,j]=2}
    if(max(c(FreqMat1[i,j],FreqMat2[i,j],FreqMat3[i,j],FreqMat4[i,j]))==FreqMat3[i,j]){FreqMat[i,j]=3}
    if(max(c(FreqMat1[i,j],FreqMat2[i,j],FreqMat3[i,j],FreqMat4[i,j]))==FreqMat4[i,j]){FreqMat[i,j]=4}
  }
}
return(FreqMat)
}
