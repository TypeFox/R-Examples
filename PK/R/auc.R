auc <- function(conc, time, group=NULL, method=c("t", "z", "boott"), alternative=c("two.sided", "less", "greater"), conf.level=0.95, strata=NULL, nsample=1000, design=c('ssd','batch','complete'), data){

  if(missing(design)) {stop("A design needs to be specified")}
  design <- match.arg(design)
  if(design=='complete'){
    if(!is.null(strata)) warning('Stratification will be ignored in a complete data design')
    auc.complete(conc=conc, time=time, group=group, method=method, alternative=alternative, conf.level=conf.level, nsample=nsample, data=data)
  }else{
    if(design=='batch'){
      if(!is.null(strata)) warning('Stratification will be ignored in a batch design')
      auc.batch(conc=conc, time=time, group=group, method=method, alternative=alternative, conf.level=conf.level, nsample=nsample, data)
    }else{
      if(!missing(conc) && ( is.list(conc) || is.list(time))) stop('Both time and concentration need to be a vector')
      auc.ssd(conc=conc, time=time, group=group, method=method, alternative=alternative, conf.level=conf.level, strata=strata, nsample=nsample, data)
    }
  }
}

