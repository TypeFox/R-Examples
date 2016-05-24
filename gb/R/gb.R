gboot <- function(x,gldobj,statistic,...)
  {
    if(!(class(gldobj)=='gld'||class(gldobj)=='egld'))
      stop("'gldobj' must be a fitted GLD/EGLD.")
    rgld <- function(data,gldobj)
      ifelse(class(gldobj)=="gld",
             rrsgld(length(data),gldobj$para),
             regld(length(data),gldobj$para))
    boot(x,statistic, sim='parametric',
         ran.gen = rgld, mle = gldobj,...)
  }
