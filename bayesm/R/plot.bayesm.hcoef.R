plot.bayesm.hcoef=function(x,names,burnin=trunc(.1*R),...){
#
# S3 method to plot arrays of draws of coefs in hier models
#   3 dimensional arrays:  unit x var x draw
#   P. Rossi 2/07
#
  X=x
  if(mode(X) == "list") stop("list entered \n Possible Fixup: extract from list \n")
  if(mode(X) !="numeric") stop("Requires numeric argument \n")
  d=dim(X)
  if(length(d) !=3) stop("Requires 3-dim array \n") 
  op=par(no.readonly=TRUE)
  on.exit(par(op))
  on.exit(devAskNewPage(FALSE),add=TRUE)
  nunits=d[1]
  nvar=d[2]
  R=d[3]
  if(missing(names)) {names=as.character(1:nvar)}
  if(R < 100) {cat("fewer than 100 draws submitted \n"); return(invisible())}
  #
  #  plot posterior distributions of nvar coef for 30 rand units
  #
 
  rsam=sort(sample(c(1:nunits),30))  # randomly sample 30 cross-sectional units
  par(mfrow=c(1,1))
  par(las=3)  # horizontal labeling
  devAskNewPage(TRUE)
  for(var in 1:nvar){
       ext=X[rsam,var,(burnin+1):R]; ext=data.frame(t(ext))
       colnames(ext)=as.character(rsam)
       out=boxplot(ext,plot=FALSE,...)
       out$stats=apply(ext,2,quantile,probs=c(0,.05,.95,1))
       bxp(out,xlab="Cross-sectional Unit",main=paste("Coefficients on Var ",names[var],sep=""),boxfill="magenta",...)
  }
  #
  # plot posterior means for each var 
  #
  par(las=1)
  pmeans=matrix(0,nrow=nunits,ncol=nvar)
  for(i in 1:nunits) pmeans[i,]=apply(X[i,,(burnin+1):R],1,mean)
  attributes(pmeans)$class="bayesm.mat"
  for(var in 1:nvar) names[var]=paste("Posterior Means of Coef ",names[var],sep="")
  plot(pmeans,names,TRACEPLOT=FALSE,INT=FALSE,DEN=FALSE,CHECK_NDRAWS=FALSE,...)
invisible()
}
