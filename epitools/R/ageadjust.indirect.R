"ageadjust.indirect" <-
  function(count, pop, stdcount, stdpop, stdrate=NULL, conf.level = 0.95){
  zv <- qnorm(0.5*(1+conf.level))
  countsum <- sum(count)
  if(is.null(stdrate)==TRUE & length(stdcount)>1 & length(stdpop>1)){
    stdrate <- stdcount/stdpop
  }
  ##indirect age standardization
  ##a. sir calculation
  expected <- sum(stdrate * pop)
  sir <- countsum/expected
  logsir.lci <- log(sir) - zv * (1/sqrt(countsum))
  logsir.uci <- log(sir) + zv * (1/sqrt(countsum))
  sir.lci <- exp(logsir.lci)
  sir.uci <- exp(logsir.uci)
                
  ##b. israte calculation
  stdcrate <- sum(stdcount)/sum(stdpop)
  crude.rate <- sum(count)/sum(pop)
  isr <- sir * stdcrate
  isr.lci <- sir.lci * stdcrate
  isr.uci <- sir.uci * stdcrate

  results <- list(sir=c(observed=countsum,exp=expected,
                    sir=sir,lci=sir.lci,uci=sir.uci), 
                  rate=c(crude.rate=crude.rate,adj.rate=isr,
                    lci=isr.lci,uci=isr.uci))
  results
}
