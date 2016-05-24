"bic.control"<-function(factor=function(nobs){log(nobs)}){
  ## Purpose: Control function for tune.ahazpen
  ##          - tuning via BIC type criterion
  ## ----------------------------------------------------------------------
  ## Arguments:
  ##   factor        : how much should df penalize? We minimize criterion 
  ##                   loss(lambda)+df(lambda)*factor(nobs)/nobs
  ## ----------------------------------------------------------------------
  ## Author: Anders Gorst-Rasmussen
  return(list("factor"=factor,"type"="BIC"))
}

"cv.control"<-function(nfolds=5,reps=1,foldid=NULL,trace=FALSE)
  {
    ## Purpose: Control function for tune.ahazpen
    ##          - tuning via cross-validation
    ## ----------------------------------------------------------------------
    ## Arguments:
    ##   nfolds        : number of folds
    ##   rep           : number of repetitions of nfolds cross-validation
    ##   foldid        : optional vector of integers describing which fold each
    ##                   observation belongs to
    ##   trace         : print progress?
    ## ----------------------------------------------------------------------
    ## Author: Anders Gorst-Rasmussen
    getfolds<-function(nobs)
      {
        if(!is.numeric(nfolds))
          stop("invalid 'nfolds'")
        if(is.null(foldid)){
          folds <- ahaz.cvfolds(nobs, nfolds)
              return(list("nfolds"=length(folds),"all.folds"=folds))
        } else {
          if(!is.numeric(foldid) || length(foldid) != nobs ||  length(unique(foldid)) != max(foldid) || sort(unique(foldid)) != 1:max(foldid))
            stop("invalid 'foldid'")
          folds <- split(1:nobs,foldid)
          nfolds <- length(folds)
          if(min(unlist(lapply(folds,length)))<=2)
            stop("too few observations in one or more folds")
          return(list("nfolds"=nfolds,"all.folds"=folds))
        }
      }
    return(list("type"="CV","getfolds"=getfolds,"trace"=trace,"rep"=reps))
  }
