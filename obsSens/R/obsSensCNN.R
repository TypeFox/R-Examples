`obsSensCNN` <-
function(model, which = 2,
                       gamma=round(seq(0,2*bstar,length=6),4), 
                       rho = c(0,0.5,0.75,0.85,0.9,0.95,0.98,0.99),
                       sdx,
                       logOdds=FALSE, method=c('approx','sim')) {

  method <- match.arg(method)
  if(method=='sim'){
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

  u <- expand.grid(rho=rho,gamma=gamma)

  if (logOdds) {
    a <- u$gamma*u$rho/sdx
    b <- bstar - a
    b.ll <- bstar.ci[1] - a
    b.ul <- bstar.ci[2] - a
  } else {
    r <- exp(bstar)
    r.ll <- exp(bstar.ci[1])
    r.ul <- exp(bstar.ci[2])
    a <- exp(u$gamma * u$rho/sdx)
    b <- r/a
    b.ll <- r.ll/a
    b.ul <- r.ul/a
  }

  out <- list()

  out$beta <- array( b, c(length(rho),length(gamma)),
                    list(rho=rho,gamma=gamma))
  out$lcl <- array( b.ll, c(length(rho),length(gamma)),
                    list(rho=rho,gamma=gamma))
  out$ucl <- array( b.ul, c(length(rho),length(gamma)),
                    list(rho=rho,gamma=gamma))
  
  out$log <- logOdds
  out$xname <- names(bstar)
  out$type <- 'cat'
  class(out) <- 'obsSens'
  
  out
}

