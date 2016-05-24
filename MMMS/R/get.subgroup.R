get.subgroup <-
function(score.info,data,lb=20,ub=80) {
  
  # get sub info for all subjects
  sub.info=score2subgroup.all(score.info$score.all,data$treat)
    
  index = which(sub.info$pct>=lb & sub.info$pct<=ub)
  if(length(index)==0) {
    return(list(pct=NA,cutoff=NA,chisq=NA,i.best=NA,fit.best=NA,
                pct.best=NA,cutoff.best=NA,chisq.best=-Inf))
  }
  
  pct=sub.info$pct[index]
  cutoff=sub.info$cutoff[index]
  
  chisq = c()
  fit0 = fit = vector("list",length(index))
  for(ii in seq_along(index)) {
    data.tmp=data
    data.tmp$sub.main=sub.info$subs[,index[ii]]
    data.tmp$sub=data.tmp$sub.main*data.tmp$treat
    data.tmp$bio.main=score.info$score.main
    
    if(length(data.tmp$covar)>0) {
      if(sd(data.tmp$bio.main)==0) {
        fit0[[ii]] = with(data.tmp,coxph(Surv(time,event)~treat+covar))
        fit[[ii]]  = with(data.tmp,coxph(Surv(time,event)~sub+treat+covar))
        chisq[ii] = 2*(fit[[ii]]$loglik-fit0[[ii]]$loglik)[2]
      } else {
        fit0[[ii]] = with(data.tmp,coxph(Surv(time,event)~treat+covar+bio.main))
        fit[[ii]]  = with(data.tmp,coxph(Surv(time,event)~sub+treat+covar+bio.main))
        chisq[ii] = 2*(fit[[ii]]$loglik-fit0[[ii]]$loglik)[2]
      }
    } else {
      if(sd(data.tmp$bio.main)==0) {
        fit0[[ii]] = with(data.tmp,coxph(Surv(time,event)~treat))
        fit[[ii]]  = with(data.tmp,coxph(Surv(time,event)~sub+treat))
        chisq[ii] = 2*(fit[[ii]]$loglik-fit0[[ii]]$loglik)[2]
      } else {
        fit0[[ii]] = with(data.tmp,coxph(Surv(time,event)~treat+bio.main))
        fit[[ii]]  = with(data.tmp,coxph(Surv(time,event)~sub+treat+bio.main))
        chisq[ii] = 2*(fit[[ii]]$loglik-fit0[[ii]]$loglik)[2]
      }
    }
  }
  
  i.best=which.max(chisq)
  
  return(list(pct=pct,cutoff=cutoff,chisq=chisq,i.best=i.best,fit.best=summary(fit[[i.best]]),
              pct.best=pct[i.best],cutoff.best=cutoff[i.best],chisq.best=chisq[i.best]))
}
