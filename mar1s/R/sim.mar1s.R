`sim.mar1s` <-
function(object, n.ahead = 1, n.sim = 1, start.time = 0,
         xreg.absdata = NULL, init.absdata = NULL)
{
  arcoef <- head(coef(object$logstoch.ar1), 1)
  xregcoef <- tail(coef(object$logstoch.ar1), -1)  
  loginnov <- matrix(rnorm(n.ahead*n.sim, sd = object$logresid.sd),
                     n.ahead, n.sim)
  d <- .decomp(object, start.time, xreg.absdata, init.absdata)

  y1 <- compose.ar1(arcoef, loginnov, head(d$init.logstoch, 1),
                    xregcoef, d$xreg.logstoch, tail(d$init.logstoch, -1))

  cycl <- cycle(ts(y1, start = start.time,
                   frequency = frequency(object$logseasonal)))
  result <- exp(tail(y1, 1) +
                as.matrix(object$logseasonal)[tail(cycl, 1), 1])
  return(as.vector(result))
}
