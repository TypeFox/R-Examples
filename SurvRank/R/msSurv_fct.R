#' Function selects for the ranking function approaches
#'
#' This wrapper function passes the ranking method to further functions.
#' @param data data input
#' @param fold outer CV fold
#' @param fs.method ranking method
#' @param ncl number of clusters for parallel computing
#' @param cv.out number of folds in outer cross validation loop (for estimation of the predictive accuracy)
#' @param cv.in number of folds in inner cross validation loop (for model selection on the training set)
#' @param nr.var Number of variables up to which stepwise selection should be carried out. Has to be smaller than n number of observations.
#' @param t current t.times argument
#' @param c.time as defined in package \code{survAUC} time; a positive number restricting the upper limit of the time range under consideration.
#' @keywords internal
#' @export

msSurv_fct = function(data,fold,fs.method,ncl,cv.out,cv.in,nr.var,t=1,sd1=0.9,c.time,...){
  if(fs.method=="lasso.rank"){
    output = fsRankLassoAUC_fct(data=data,fold=fold,ncl=ncl,cv.out=cv.out,cv.in=cv.in,nr.var=nr.var,t=t,c.time=c.time,...)
  }
  if(fs.method=="conc.rank"){
    output = fsRankConcAUC_fct(data=data,fold=fold,ncl=ncl,cv.out=cv.out,cv.in=cv.in,nr.var=nr.var,t=t,c.time=c.time,...)
  }
  if(fs.method=="rf.rank"){
    output = fsRankRfAUC_fct(data=data,fold=fold,ncl=ncl,cv.out=cv.out,cv.in=cv.in,nr.var=nr.var,t=t,c.time=c.time,...)
  }
  if(fs.method=="cox.rank"){
    output = fsRankCoxAUC_fct(data=data,fold=fold,ncl=ncl,cv.out=cv.out,cv.in=cv.in,nr.var=nr.var,t=t,c.time=c.time,...)
  }
  if(fs.method=="rf.rank"){
    output = fsRankRfAUC_fct(data=data,fold=fold,ncl=ncl,cv.out=cv.out,cv.in=cv.in,nr.var=nr.var,t=t,c.time=c.time,...)
  }
  if(fs.method=="boost.rank"){
    output = fsRankBoostAUC_fct(data=data,fold=fold,ncl=ncl,cv.out=cv.out,cv.in=cv.in,nr.var=nr.var,t=t,c.time=c.time,...)
  }
  if(fs.method=="rpart.rank"){
    output = fsRankRpartAUC_fct(data=data,fold=fold,ncl=ncl,cv.out=cv.out,cv.in=cv.in,nr.var=nr.var,t=t,c.time=c.time,...)
  }
  if(fs.method=="randcox.rank"){
    output = fsRankMonsterAUC_fct(data=data,fold=fold,ncl=ncl,cv.out=cv.out,cv.in=cv.in,nr.var=nr.var,t=t,c.time=c.time,...)
  }
  if(fs.method=="wang.rank"){
    output = fsRankWangAUC_fct(data=data,fold=fold,ncl=ncl,cv.out=cv.out,cv.in=cv.in,nr.var=nr.var,t=t,c.time=c.time,...)
  }
  return(output)
}
