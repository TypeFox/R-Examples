get.score <-
function(time,event,treat,bio,covar=NULL,nfolds=5,alpha=0.5,pos.direction=FALSE) {
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
  x=cbind(treat,covar,bio,c(treat)*bio) 
  index.treat=1 # treat
  index.bio  =1+n.covar+(1:n.bio) # biomarkers
  index.inter=1+n.covar+n.bio+(1:n.bio) # treat*biomarkers
  pf = rep(0,n.vars)
  pf[c(index.bio,index.inter)]=1
  
  # glmnet regression
  lam = cv.glmnet(x=x,y=surv,family="cox",alpha=alpha,standardize=F,
                  penalty.factor=pf,nlambda=200,nfolds=nfolds)$lambda.min
  fit = glmnet(x=x,y=surv,family="cox",alpha=alpha,standardize=F,penalty.factor=pf,nlambda=200)
  
  lam.best= fit$lambda[which.min(abs(fit$lambda-lam))]
  
  nvars.min=1
  coefs = coef(fit,s=lam.best)[index.inter,1]
  coefs.main = coef(fit,s=lam.best)[index.bio,1]
  
  betas=as(fit$beta,'matrix')
  
  if(sum(coefs!=0)<nvars.min) {
    n.diff = colSums(betas[index.inter,,drop=FALSE]!=0)-nvars.min # ideally 0, or smallest positive
    n.diff[n.diff<0] = NA
    
    if(sum(!is.na(n.diff))==0) {
      sma.p=c()
      sma.coef=c()
      for(i in 1:n.bio) {
        if(sd(bio[,i])==0 | sd(bio[treat==1,i])==0) {
          sma.p[i]=1
          sma.coef[i]=0
        } else {
          if(n.covar>0) {
            fit=summary(coxph(surv~I(treat*bio[,i])+treat+covar+bio[,i]))
          } else {
            fit=summary(coxph(surv~I(treat*bio[,i])+treat+bio[,i]))
          }
          sma.p[i]=fit$coef['I(treat * bio[, i])','Pr(>|z|)']
          sma.coef[i]=fit$coef['I(treat * bio[, i])','coef']
        }
      }
      coefs[which.min(sma.p)]=sma.coef[which.min(sma.p)]
      lam.best=NA
    } else {
      coefs = betas[index.inter,which.min(n.diff)]
      lam.best = fit$lambda[which.min(n.diff)]
    }
  }
  
  score.all = ifelse(pos.direction,1,-1) * c(bio %*% coefs)
  score = score.all[treat==1]
  
  score.main = c(bio %*% coefs.main)
  
  return(list(score=score,score.all=score.all,score.main=score.main,
              coefs=coefs,coefs.main=coefs.main,fit=fit,lam.best=lam.best,
              treat=treat,alpha=alpha))
}
