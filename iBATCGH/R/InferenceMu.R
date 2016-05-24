InferenceMu <-
function(listMu,niter,burnin){
mumean=rep(0,4)

for(i in burnin:niter){
  temp=listMu[[i]]
  mumean=mumean+temp
}
mumean=mumean/(niter-burnin)
return(mumean)
}
