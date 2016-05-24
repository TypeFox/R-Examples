summary.bayesm.nmix=function(object,names,burnin=trunc(.1*nrow(probdraw)),...){
  nmixlist=object
  if(mode(nmixlist) != "list") stop(" Argument must be a list \n")
  probdraw=nmixlist[[1]]; compdraw=nmixlist[[3]]
  if(!is.matrix(probdraw)) stop(" First Element of List (probdraw) must be a matrix \n")
  if(mode(compdraw) != "list") stop(" Third Element of List (compdraw) must be a list \n")
  ncomp=length(compdraw[[1]])
  if(ncol(probdraw) != ncomp) stop(" Dim of First Element of List not compatible with Dim of Second
      \n")
#
# function to summarize draws of normal mixture components
#
R=nrow(probdraw)
if(R < 100) {cat("fewer than 100 draws submitted \n"); return(invisible())}
datad=length(compdraw[[1]][[1]]$mu)
mumat=matrix(0,nrow=R,ncol=datad)
sigmat=matrix(0,nrow=R,ncol=(datad*datad))
if(missing(names)) names=as.character(1:datad)
for(i in (burnin+1):R){
   if(i%%500 ==0) cat("processing draw ",i,"\n",sep="");fsh()
   out=momMix(probdraw[i,,drop=FALSE],compdraw[i])
   mumat[i,]=out$mu
   sigmat[i,]=out$sigma
   } 
cat("\nNormal Mixture Moments\n Mean\n")
attributes(mumat)$class="bayesm.mat"
attributes(sigmat)$class="bayesm.var"
summary(mumat,names,burnin=burnin,QUANTILES=FALSE,TRAILER=FALSE)
cat(" \n")
summary(sigmat,burnin=burnin)
cat("note: 1st and 2nd Moments for a Normal Mixture \n")
cat("      may not be interpretable, consider plots\n")
invisible()
}

