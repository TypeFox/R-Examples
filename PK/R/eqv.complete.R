eqv.complete <- function(conc, time, group, dependent=FALSE, method=c("fieller", "z", "boott"), conf.level=0.90, nsample=1000, data){
  message('Confidence intervals for AUCs are found using a batch design with one batch. Please see the help file "auc.complete" for some alternative examples.')

  if(!missing(data)){
    cnames <- colnames(data)
    if(!any(cnames=='conc')){stop("data does not contain a variable conc")}
    if(!any(cnames=='time')){stop("data does not contain a variable time")}
    conc <- data$conc
    time <- data$time
    if(any(cnames=='group')){
      group <- data$group
    }
  }

  time <- list(batch1=unlist(time))
  conc <- list(batch1=unlist(conc))
  if(!is.null(group)){
    group <- list(batch1=unlist(group))
  }

  res <- eqv.batch(conc=conc, time=time, group=group, dependent=dependent, method=method, conf.level=conf.level, nsample=nsample)
  res$design <- "complete"
  return(res)
}
