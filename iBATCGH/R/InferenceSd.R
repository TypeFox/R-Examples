InferenceSd <-
function(listSd,niter,burnin){
sigmamean=rep(0,4)

for(i in burnin:niter){
  temp=listSd[[i]]
  sigmamean=sigmamean+temp
}
sigmamean=sigmamean/(niter-burnin)
return(sigmamean)
}
