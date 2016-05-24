get.score.main <-
function(time,event,treat,bio,covar=NULL,nfolds=5,alpha=0.5) {
  # check input data
  stopifnot(length(time)==length(event))
  stopifnot(length(time)==length(treat))
  stopifnot(!missing(bio))  
  stopifnot(treat %in% c(0,1))
  stopifnot(sd(c(treat))>0)
  stopifnot(event %in% c(0,1))
  
  if(any(is.na(time))) stop('No missing value is expected in time.')
  if(any(is.na(event))) stop('No missing value is expected in event.')
  if(any(is.na(treat))) stop('No missing value is expected in treat.')
  if(any(is.na(bio))) stop('No missing value is expected in bio.')
  if(!is.null(covar)) { 
    if(any(is.na(covar))) stop('No missing value is expected in covar.') 
    if(!is.matrix(covar)) stop('covar needs to be a numeric matrix.') 
  }
  
  if(is.null(covar)) {
    n.covar=0 
  } else {
    covar=cbind(covar)
    n.covar=ncol(covar)
    stopifnot(n.covar>0)
    stopifnot(nrow(covar)==length(time))
  }
  
  bio=cbind(bio)
  n.bio=ncol(bio)
  
  # numbers
  n.vars=1+n.covar+n.bio+n.bio
  n=length(treat)
  
  # survival data
  surv = Surv(time,event)
  
  # construct x
  x=cbind(treat,covar,bio) 
  index.treat=1 # treat
  index.bio  =1+n.covar+(1:n.bio) # biomarkers
  pf = rep(0,n.vars)
  pf[index.bio]=1
  
  # glmnet regression
  lam = cv.glmnet(x=x,y=surv,family="cox",alpha=alpha,standardize=FALSE,
                  penalty.factor=pf,nlambda=200,nfolds=nfolds)$lambda.min
  fit = glmnet(x=x,y=surv,family="cox",alpha=alpha,standardize=F,penalty.factor=pf,nlambda=200)
  
  lam.best= fit$lambda[which.min(abs(fit$lambda-lam))]
  
  coefs = coef(fit,s=lam.best)[,1]
  index.selected = abs(coefs)>0
  coefs.selected = coefs[index.selected]
  
  xx = x[,index.selected]
  fit.selected = coxph(surv~xx,init=coefs.selected,iter=0)
  sfit = survfit(fit.selected,newdata=as.data.frame(xx))
  
  return(list(fit=fit,lam.best=lam.best,fit.selected=fit.selected,sfit=sfit))
}
