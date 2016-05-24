# boot : boot

bootstrap <-
function(x,fun,nrep=1000,conf.level=0.95,...){
  simul <- boot::boot(x,fun,R=nrep,...)
  estimate <- simul$t0
  names(estimate) <- "original value"
  interval <- .ci(simul$t,conf.level=conf.level)
  attr(interval,"conf.level") <- conf.level
  dname <- paste(deparse(substitute(x)),"\n",nrep," replicates",sep="")
  result <- list(method="Bootstrap",data.name=dname,estimate=estimate,conf.level=conf.level,rep=nrep,conf.int=interval)
  class(result) <- "htest"
  return(result)
}
