infExpoTree <- function(pars,times,ttypes,survival=TRUE,vflag=0,
                        root.lineages=0) 
{
  f <- .Call("infTreeEval",parameters=pars,
      times=as.numeric(times),ttypes=as.integer(ttypes),
      survival=as.integer(c(survival,vflag,root.lineages)))
  return(f)
}

