`obsSensSCN` <-
function(model, which = 2,
                       gamma=round(seq(0,3,length=6),4),
                       delta = seq(0,3,0.5),
                       logHaz=FALSE, method=c('approx','sim')) {
  
  method <- match.arg(method)
  if(method=='sim') {
    stop('simulation method not implemented yet')
  }

  out <- obsSensCCN(model=model, which=which, gamma=gamma,
               delta=delta, logOdds=logHaz, method=method)
  
  out$type <- 'surv'
  
  out
}

