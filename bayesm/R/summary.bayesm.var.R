summary.bayesm.var=function(object,names,burnin=trunc(.1*nrow(Vard)),tvalues,QUANTILES=FALSE,...){
#
# S3 method to summarize draws of var-cov matrix (stored as a vector)
#   Vard is R x d**2 array of draws
#   P. Rossi 2/07
#
  Vard=object
  if(mode(Vard) == "list") stop("list entered \n Possible Fixup: extract from list \n")
  if(!is.matrix(Vard)) stop("Requires matrix argument \n")
  if(trunc(sqrt(ncol(Vard)))!=sqrt(ncol(Vard))) stop("Argument cannot be draws from a square matrix \n")
  if(nrow(Vard) < 100) {cat("fewer than 100 draws submitted \n"); return(invisible())}
  d=sqrt(ncol(Vard))
  corrd=t(apply(Vard[(burnin+1):nrow(Vard),],1,nmat))
  pmeancorr=apply(corrd,2,mean)
  dim(pmeancorr)=c(d,d)
  indexdiag=(0:(d-1))*d+1:d
  var=Vard[(burnin+1):nrow(Vard),indexdiag]
  sdd=sqrt(var)
  pmeansd=apply(sdd,2,mean)
  mat=cbind(pmeansd,pmeancorr)
  if(missing(names)) names=as.character(1:d)

  cat("Posterior Means of Std Deviations and Correlation Matrix \n")
  rownames(mat)=names
  colnames(mat)=c("Std Dev",names)
  print(mat,digits=2)

  cat("\nUpper Triangle of Var-Cov Matrix \n")
  ind=as.vector(upper.tri(matrix(0,ncol=d,nrow=d),diag=TRUE))
  labels=cbind(rep(c(1:d),d),rep(c(1:d),each=d))
  labels=labels[ind,]
  plabels=paste(labels[,1],labels[,2],sep=",")
  uppertri=as.matrix(Vard[,ind])
  attributes(uppertri)$class="bayesm.mat"
  summary(uppertri,names=plabels,burnin=burnin,tvalues=tvalues,QUANTILES=QUANTILES)
  invisible()
}
