predValues <- function(sens, spec, prev){
  if(all(sens < 0 | sens > 1))
    stop("'sens' has to be in [0,1]")
  if(all(spec < 0 | spec > 1))
    stop("'spec' has to be in [0,1]")
  if(length(sens) != length(spec))
    stop("lengths of 'sens' and 'spec' have to be identical")
  if(length(prev) > 1)
    stop("'pre' has to be of length 1")
  if(prev < 0 | prev > 1)
    stop("'pre' has to be in [0,1]")
  
  ppv <- sens*prev/(sens*prev + (1-spec)*(1-prev))
  npv <- spec*(1-prev)/(spec*(1-prev) + (1-sens)*prev)
  if(length(sens) > 1)
    res <- cbind("PPV" = ppv, "NPV" = npv)
  else
    res <- c("PPV" = ppv, "NPV" = npv)
  res
}
