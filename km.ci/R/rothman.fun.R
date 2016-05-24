"rothman.fun" <-
function(sfit,conf.level=0.95)
{
    if(conf.level < 0 || conf.level > 1)
      stop("confidence level must be between 0 and 1")
    else
    s.t <- sfit$surv
    zalpha<-qnorm( 1 - (1-conf.level)/2)
    tvec<-sfit$n.event/(sfit$n.risk*(sfit$n.risk - sfit$n.event))
    tv2<-cumsum(tvec)
    var.st <- tv2*s.t^2
    n.null <- s.t*(1-s.t)/var.st
    roth.upper <- n.null/(n.null+zalpha^2)*(s.t+zalpha^2/(2*n.null)+zalpha*sqrt(var.st+zalpha^2/(4*n.null^2)))
    roth.lower <- n.null/(n.null+zalpha^2)*(s.t+zalpha^2/(2*n.null)-zalpha*sqrt(var.st+zalpha^2/(4*n.null^2)))
    sfit <- sfit
    roth.upper[is.na(roth.upper)] <- 1
    roth.lower[is.na(roth.lower)] <- 1
    sfit$upper <- roth.upper
    sfit$lower <- roth.lower
    return(list(rothman.upper=roth.upper,rothman.lower=roth.lower,surv.object=sfit))
}

