"ageadjust.direct" <-
  function (count, pop, rate = NULL, stdpop, conf.level = 0.95) 
{
  if (missing(count) == TRUE & !missing(pop) == TRUE & is.null(rate) == 
      TRUE) 
    count <- rate * pop
  if (missing(pop) == TRUE & !missing(count) == TRUE & is.null(rate) == 
      TRUE) 
    pop <- count/rate
  if (is.null(rate) == TRUE & !missing(count) == TRUE & !missing(pop) == 
      TRUE) 
    rate <- count/pop
  alpha <- 1 - conf.level
  cruderate <- sum(count)/sum(pop)
  stdwt <- stdpop/sum(stdpop)
  dsr <- sum(stdwt * rate)
  dsr.var <- sum((stdwt^2) * (count/pop^2))
  wm<- max(stdwt/pop)
  gamma.lci <- qgamma(alpha/2, shape = (dsr^2)/dsr.var, scale = dsr.var/dsr)
  gamma.uci <- qgamma(1 - alpha/2, shape = ((dsr+wm)^2)/(dsr.var+wm^2), 
                      scale = (dsr.var+wm^2)/(dsr+wm))
  c(crude.rate = cruderate, adj.rate = dsr, lci = gamma.lci, 
    uci = gamma.uci)
}
