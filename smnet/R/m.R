m<-function (...,  k = -1) 
{
  cyclic <- F
  vars   <- as.list(substitute(list(...)))[-1]
  d      <- length(vars)
  term   <- deparse(vars[[1]], backtick = TRUE, width.cutoff = 500)
  
  
  if (d > 1) for (i in 2:d) term[i] <- deparse(vars[[i]], backtick = TRUE, width.cutoff = 500)
  for (i in 1:d)            term[i] <- attr(terms(reformulate(term[i])), "term.labels")
  
  # some warning messages 
  if (all.equal(round(k, 0), k) != TRUE) stop("argument k of m() should be integer")
  if (length(unique(term)) != d)         stop("Repeated variables as arguments of a smooth are not permitted")
  
  
  full.call <- paste("m(", term[1], sep = "")
  if (d > 1) for (i in 2:d) full.call <- paste(full.call, ",", term[i], sep = "")
  label <- paste(full.call, ")", sep = "")
  list(term = term, bs.dim = k,  cyclic = cyclic)
}