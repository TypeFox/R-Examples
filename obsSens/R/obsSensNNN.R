`obsSensNNN` <-
function(model, which = 2,
                       gamma=round(seq(0,2*bstar,length=6),4),
                       rho = c(0,0.5,0.75,0.85,0.9,0.95,0.98,0.99), sdx,
                       log=TRUE, method=c('approx','sim')) {
  method <- match.arg(method)
  if(method=='sim') {
    stop('simulation method not implemented yet')
  }

 if(is.numeric(model)){
    bstar <- model
    bstar.ci <- which
  } else {
    bstar <- coef(model)[which]
    bstar.ci <- confint(model, which)
    sdx <- sd( model.matrix(model)[,which] )
  }

  out <- obsSensCNN(model=model, which=which, gamma=gamma,
               rho=rho, sdx=sdx, logOdds=log, method=method)

  out$type <- 'num'

  out
}

