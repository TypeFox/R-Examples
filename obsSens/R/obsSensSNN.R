`obsSensSNN` <-
function(model, which = 2,
                       gamma=round(seq(0,3,length=6),4),
                       rho = c(0,0.5,0.75,0.85,0.9,0.95,0.98,0.99), sdx,
                       logHaz=FALSE, method=c('approx','sim')) {

  method <- match.arg(method)
  if(method=='sim') {
    stop('simulation method not implemented yet')
  }

  out <- obsSensCNN(model=model, which=which, gamma=gamma,
               rho=rho, sdx=sdx, logOdds=logHaz, method=method)

  out$type <- 'surv'

  out
}

