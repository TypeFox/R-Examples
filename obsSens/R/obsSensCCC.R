`obsSensCCC` <-
function(model, which = 2, g0 = c(2,6,10), 
                       g1, p0= seq(0,1,.2),
                       p1=p0, logOdds=FALSE, method=c('approx','sim')) {
  
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
        
  if(missing(g1)){ 
    u <- expand.grid( p1=p1, p0=p0, g0=g0 )
    u$g1 <- u$g0
  } else {
    u <- expand.grid( p1=p1, p0=p0, g0=g0, g1=g1 )
  }

  if (logOdds) {
    a <- log( ( u$g1 * u$p1 + (1-u$p1))/ (u$g0 * u$p0 + (1-u$p0)) )
    b <- bstar - a
    b.ll <- bstar.ci[1] - a
    b.ul <- bstar.ci[2] - a
  } else {
    r <- exp(bstar)
    r.ll <- exp(bstar.ci[1])
    r.ul <- exp(bstar.ci[2])
    a <- (u$g1 * u$p1 + (1-u$p1))/ (u$g0 * u$p0 + (1-u$p0))
    b <- r/a
    b.ll <- r.ll/a
    b.ul <- r.ul/a
  }

  out <- list()
        
  if( missing(g1) ){
    out$beta <- array( b, c(length(p1),length(p0),length(g0)),
                      list(P1=p1,P0=p0,Gamma=g0) )
    out$lcl <- array( b.ll, c(length(p1),length(p0),length(g0)),
                     list(P1=p1,P0=p0,Gamma=g0) )
    out$ucl <- array( b.ul, c(length(p1),length(p0),length(g0)),
                     list(P1=p1,P0=p0,Gamma=g0) )
  } else {
    out$beta <- array( b, c(length(p1),length(p0),length(g0),length(g1)),
                      list(P1=p1,P0=p0,Gamma0=g0,Gamma1=g1) )
    out$lcl <- array(b.ll,c(length(p1),length(p0),length(g0),length(g1)),
                     list(P1=p1,P0=p0,Gamma0=g0,Gamma1=g1) )
    out$ucl <- array(b.ul,c(length(p1),length(p0),length(g0),length(g1)),
                     list(P1=p1,P0=p0,Gamma0=g0,Gamma1=g1) )
  }

  out$log <- logOdds
  out$xname <- names(bstar)
  out$type <- 'cat'
  class(out) <- 'obsSens'

  out
}

