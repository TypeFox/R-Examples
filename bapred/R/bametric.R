bametric <-
function(xba, batch, y, x, metric=c("sep", "avedist", 
  "kldist", "skew", "pvca", "diffexpr", "cor"), method, ...) {
  
if(!(metric %in% c("sep", "avedist", 
  "kldist", "skew", "pvca", "diffexpr", "cor")))
    stop("Input parameter 'metric' has to be one of the following:\n'sep', 'avedist', 'kldist', 'skew', 'pvca', 'diffexpr', 'cor'.")
    
if(metric=="sep")
  return(sepscore(xba=xba, batch=batch, ...))

if(metric=="avedist")
  return(avedist(xba=xba, batch=batch))

if(metric=="kldist")
  return(kldist(xba=xba, batch=batch))

if(metric=="skew")
  return(skewdiv(xba=xba, batch=batch))

if(metric=="pvca")
  return(pvcam(xba=xba, batch=batch, y=y, ...))

if(metric=="diffexpr")
  return(diffexprm(x=x, batch=batch, y=y, method=method))

if(metric=="cor")
  return(corba(xba=xba, x=x))

}
