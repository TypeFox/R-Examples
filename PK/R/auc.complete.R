auc.complete <- function(conc, time, group=NULL, method=c("t", "z", "boott"),  alternative=c("two.sided", "less", "greater"), conf.level=0.95, nsample=1000, data) {

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

  if(length(unique(time))<length(time)) {
  message('Confidence intervals for AUCs are found using a batch design with one batch. Please see the help file "auc.complete" for some alternative examples.')
  }


  time <- list(batch1=unlist(time))
  conc <- list(batch1=unlist(conc))
  if(!is.null(group)){
    group <- list(batch1=unlist(group))
  }

  ## check if at least 2 subjects are provided. If not make sure that bootstrap method is not used.
  if((table(time)[1]==1 & is.null(group)) || (table(time)[1]==2 & length(unique(group))==2) ) {
    method <- "t"
  }
  res <- auc.batch(conc=conc, time=time, group=group, method=method, alternative=alternative, conf.level=conf.level, nsample=nsample)
  if(any(res$CIs["stderr"]==0)){
    res$CIs[c("stderr","lower","upper","df")] <- NA
  }
  res$design <- "complete"
  return(res)
}
