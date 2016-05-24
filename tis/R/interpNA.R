interpNA <- function(x, method = "constant", useTimes = F, offset = 1,
                     rule = 2, f = 0, ...){
  ## interpolates missing values in columns of tis x
  if(!inherits(x, "tis"))
    stop(paste(deparse(substitute(x)), "is not a tis object"))
  xtimes <- xti <- ti(x)
  if(useTimes) xtimes <- time(xti, offset = offset)
  
  ## internal function interpVec does all the real work
  interpVec <- function(xc){  ## xc is a column of x
    naSpots <- is.na(xc)
    if(length(naSpots) > 0){
      naTimes <- xtimes[naSpots]
      interpfun <- switch(method,
                          constant = , linear = {
                            approxfun(xtimes, y = xc, method = method,
                                      ..., f = f, rule = rule)
                          },
                          fmm = , natural = , periodic = {
                            splinefun(xtimes, y = xc, method = method)
                          },
                          stop(paste("unknown method", method)))
      xc[naSpots] <- interpfun(naTimes)
    }
    xc
  }
  
  if(is.matrix(x))
    z <- apply(unclass(x), 2, interpVec)
  else
    z <- interpVec(x)
  attributes(z) <- attributes(x)
  z
}
