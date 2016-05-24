#' Wrapper function using the univariate Cox models as ranking method
#'
#' This wrapper function passes the ranking method to further functions. Setting up for parallel computing of different folds via foreach.
#' @param input see details in \code{\link{CVrankSurv_fct}}
#' @param ... other arguments, not used now
#' @keywords internal
#' @export

fsRankCoxAUC_fct = function(data,fold,ncl,cv.out,cv.in,nr.var,t=1,sd1=0.9,c.time,...){
  data.out = ranking = pred.in = pred.out = pred.out1 = used.rank = used.rank1se = list()
  auc.out = auc.out1 = NULL
  out.s=matrix(NA,nrow=cv.in,ncol=nr.var)
  i=NULL
  #####
  cl=parallel::makeCluster(ncl,type="PSOCK")
  doParallel::registerDoParallel(cl)
  output = foreach::foreach(i = 1:cv.out, .packages=c('SurvRank'),.errorhandling = "pass",.verbose = F) %dopar% {
    innerAUC_fct(f=fsSurvRankCox,data=data,t=t,cv.out=cv.out,cv.in=cv.in,i=i,fold=fold,data.out=data.out,nr.var=nr.var,out.s=out.s,sd1=sd1,c.time=c.time,
                 ranking=ranking,used.rank=used.rank,used.rank1se=used.rank1se,pred.in=pred.in,pred.out=pred.out,pred.out1=pred.out1,auc.out=auc.out,auc.out1=auc.out1)
  }
  parallel::stopCluster(cl)
  res=output
  return(res)
}
