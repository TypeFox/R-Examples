MMMS <-
function(time,event,treat,bio,covar=NULL,pct.lb=20,pct.ub=80,n.boot=1000,
               pos.direction=FALSE,nfolds=5,alpha=0.5,verbose=FALSE,seed=NULL) {
  if(any(is.na(time))) stop('No missing value is expected in time.')
  if(any(is.na(event))) stop('No missing value is expected in event.')
  if(any(is.na(treat))) stop('No missing value is expected in treat.')
  if(any(is.na(bio))) stop('No missing value is expected in bio.')
  if(!is.null(covar)) { 
    if(any(is.na(covar))) stop('No missing value is expected in covar.') 
    if(!is.matrix(covar)) stop('covar needs to be a numeric matrix.') 
  }
  
  # set seed if provided
  if(!is.null(seed)) set.seed(seed)
  
  # get composite scores
  score = get.score(time,event,treat,bio,covar,nfolds=nfolds,alpha=alpha,
                    pos.direction=pos.direction)
  
  data = list(time=time,event=event,treat=treat,bio=bio,covar=covar)
  
  # get subgroup
  subgrp = get.subgroup(score,data,lb=pct.lb,ub=pct.ub)
  
  # get bootstrap p-value
  p.boot = NA
  
  if(n.boot>0) {
    # set seed if provided
    if(!is.null(seed)) set.seed(seed)
    
    score.boot = vector("list",n.boot)
    sub.boot = vector("list",n.boot)
    
    main.only=get.score.main(time,event,treat,bio,covar,nfolds=nfolds,alpha=alpha)
    
    for(i in 1:n.boot) {
      simboot = sample.tte(main.only$sfit$time,main.only$sfit$surv)
      simdat.boot = data
      simdat.boot$time = simboot$times
      simdat.boot$event = simboot$event
      score.boot[[i]] = with(simdat.boot, get.score(time,event,treat,bio,covar,nfolds=nfolds,
                                                    alpha=alpha,pos.direction=pos.direction))
      sub.boot[[i]] = get.subgroup(score.boot[[i]],simdat.boot,lb=pct.lb,ub=pct.ub)
      
      if(verbose) print(i)
    }
    p.boot=get.p.boot(subgrp$chisq.best,sub.boot)  
  }
 
  return(list(score.obj=score,
              score=score$score.all,
              score.main=score$score.main,
              coefs=score$coefs,
              coefs.main=score$coefs.main,
              fit=score$fit,
              lambda=score$lam.best,
              alpha=alpha,
              subgrp.obj=subgrp,
              subgrp.size=subgrp$pct.best,
              subgrp.fit=subgrp$fit.best,
              subgrp.cut=subgrp$cutoff.best,
              subgrp.pval=p.boot,
              n.boot=n.boot))
}
