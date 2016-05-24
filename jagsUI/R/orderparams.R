
order.params <- function(samples,parameters.to.save,DIC,verbose=TRUE){
  
  params <- colnames(samples[[1]])
  params <- params[order(match(sapply(strsplit(params, "\\["), "[", 1),
                               sapply(strsplit(parameters.to.save, "\\["), "[", 1)))]
  
  if(DIC&&('deviance'%in%params)){
    params <- c(params[params!='deviance'],'deviance')
  } else if (DIC&&!('deviance'%in%params)){
    if(verbose){warning('JAGS did not monitor deviance.')}
    DIC <- FALSE
  } 
  
  samples <- samples[,params]
  
  return(samples)
  
}