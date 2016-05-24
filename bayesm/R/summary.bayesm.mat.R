summary.bayesm.mat=function(object,names,burnin=trunc(.1*nrow(X)),tvalues,QUANTILES=TRUE,TRAILER=TRUE,...){
#
# S3 method to compute and print posterior summaries for a matrix of draws
#  P. Rossi 2/07
#
  X=object
  if(mode(X) == "list") stop("list entered \n Possible Fixup: extract from list \n")
  if(mode(X) !="numeric") stop("Requires numeric argument \n")
  if(is.null(attributes(X)$dim)) X=as.matrix(X)
  nx=ncol(X)
  if(missing(names)) names=as.character(1:nx)
  if(nrow(X) < 100) {cat("fewer than 100 draws submitted \n"); return(invisible())}
  X=X[(burnin+1):nrow(X),,drop=FALSE]
  mat=matrix(apply(X,2,mean),nrow=1)
  mat=rbind(mat,sqrt(matrix(apply(X,2,var),nrow=1)))
  num_se=double(nx); rel_eff=double(nx); eff_s_size=double(nx)
  for(i in 1:nx) 
     {out=numEff(X[,i])
      if(is.nan(out$stderr)) 
          {num_se[i]=-9999; rel_eff[i]=-9999; eff_s_size[i]=-9999} 
      else
          {num_se[i]=out$stderr; rel_eff[i]=out$f; eff_s_size[i]=nrow(X)/ceiling(out$f)}
     }
  mat=rbind(mat,num_se,rel_eff,eff_s_size)
  colnames(mat)=names
  rownames(mat)[1]="mean"
  rownames(mat)[2]="std dev"
  rownames(mat)[3]="num se"
  rownames(mat)[4]="rel eff"
  rownames(mat)[5]="sam size"
  if(!missing(tvalues))
    {if(mode(tvalues)!="numeric") stop("true values arguments must be numeric \n")
     if(length(tvalues) != nx) stop("true values argument is wrong length \n")
     mat=rbind(tvalues,mat) }
  cat("Summary of Posterior Marginal Distributions ")
  cat("\nMoments \n")
  print(t(mat),digits=2)
  if(QUANTILES){
     qmat=apply(X,2,quantile,probs=c(.025,.05,.5,.95,.975))
     colnames(qmat)=names
     if(!missing(tvalues))
       { qmat=rbind(tvalues,qmat)}
     cat("\nQuantiles \n")
     print(t(qmat),digits=2)}
  if(TRAILER) cat(paste("   based on ",nrow(X)," valid draws (burn-in=",burnin,") \n",sep=""))
  invisible(t(mat))
}

