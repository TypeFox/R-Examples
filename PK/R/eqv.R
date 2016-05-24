eqv <- function(conc, time, group, dependent=FALSE, method=c("fieller", "z", "boott"), conf.level=0.90, strata=NULL, nsample=1000, design=c("ssd","batch","complete"), data) {
  if(missing(design)) {stop("A design needs to be specified")}
  design <- match.arg(design)
  if(design=='complete'){
    if(!is.null(strata)) warning('Stratification will be ignored in a complete data design')
    eqv.complete(conc=conc, time=time, group=group, dependent=dependent, method=method, conf.level=conf.level, nsample=nsample, data)
  }else{
    if(design=='batch'){
      if(!is.null(strata)) warning('Stratification will be ignored in a batch design')
      eqv.batch(conc=conc, time=time, group=group, dependent=dependent, method=method, conf.level=conf.level, nsample=nsample, data)
    }else{
      if(!missing(conc) && (is.list(conc) || is.list(time))) stop('Both time and concentration need to be a vector')
      eqv.ssd(conc=conc, time=time, group=group, dependent=dependent, method=method, conf.level=conf.level, strata=strata, nsample=nsample, data)
    }
  }
}
