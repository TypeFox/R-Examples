`obsSensSCC` <-
function(model, which = 2, g0 = c(2,6,10),
                       g1, p0 = seq(0,1,.2), p1=p0,
                       logHaz=FALSE, method=c('approx','sim')) {
  method <- match.arg(method)
  if(method=='sim') {
    stop('simulation method not implemented yet')
  }

  out <- if(missing(g1)){
    obsSensCCC(model=model, which=which, g0=g0,
               p0=p0, p1=p1, logOdds=logHaz, method=method)
  } else {
    obsSensCCC(model=model, which=which, g0=g0, g1=g1,
               p0=p0, p1=p1, logOdds=logHaz, method=method)
  }

  out$type <- 'surv'
  
  out
}

