InferenceA <-
function(listA,niter,burnin){
Amean=matrix(0,nrow=4,ncol=4)

for(i in burnin:niter){
  temp=listA[[i]]
  Amean=Amean+temp
}
Amean=Amean/(niter-burnin)
return(Amean)
}
