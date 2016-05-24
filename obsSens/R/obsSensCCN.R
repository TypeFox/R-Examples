`obsSensCCN` <-
function(model, which = 2,
                       gamma=round(seq(0,2*bstar,length=6),4), 
                       delta = seq(0,3,0.5),
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
  }

  u <- expand.grid(delta=delta,gamma=gamma)

  if (logOdds) {
    a <- u$gamma*u$delta
    b <- bstar - a
    b.ll <- bstar.ci[1] - a
    b.ul <- bstar.ci[2] - a
  } else {
    r <- exp(bstar)
    r.ll <- exp(bstar.ci[1])
    r.ul <- exp(bstar.ci[2])
    a <- exp(u$gamma * u$delta)
    b <- r/a
    b.ll <- r.ll/a
    b.ul <- r.ul/a
  }

  out <- list()

  out$beta <- array( b, c(length(delta),length(gamma)),
                    list(delta=delta,gamma=gamma))
  out$lcl <- array( b.ll, c(length(delta),length(gamma)),
                    list(delta=delta,gamma=gamma))
  out$ucl <- array( b.ul, c(length(delta),length(gamma)),
                    list(delta=delta,gamma=gamma))
  
  out$log <- logOdds
  out$xname <- names(bstar)
  out$type <- 'cat'
  class(out) <- 'obsSens'
  
  out
}

