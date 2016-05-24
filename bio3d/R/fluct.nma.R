"fluct.nma" <-
  function(nma, mode.inds=NULL) {
    
    kb <- 0.00831447086363271
    pi <- 3.14159265359
    
    if(!"nma" %in% class(nma))
      stop("fluct.nma: must supply 'nma' object, i.e. from 'nma'")
    
    if("VibrationalModes" %in% class(nma))
      mass <- TRUE
    else
      mass <- FALSE

    if(is.null(mode.inds))
      mode.inds <- seq(nma$triv.modes+1, length(nma$L))

    if(min(mode.inds)<=nma$triv.modes)
      stop("'mode.inds' should not contain indices to trivial modes")

    f <- apply(nma$U, 2, function(x) { rowSums(matrix(x, ncol=3, byrow=TRUE)**2) })
    
    if(mass)
      freq <- nma$frequencies**2
    else
      freq <- nma$force.constants


    for ( i in mode.inds )  {
      f[,i] <- f[,i] / freq[i]
    }

    if(length(mode.inds)>1)
      f <- rowSums(f[,mode.inds])
    else
      f <- f[,mode.inds]
  
    if(mass) {
      f <- f / nma$mass
      s <- 1/(2*pi)**2
      if(!is.null(nma$temp))
        s <- s*kb*nma$temp
      f <- f*s
    }
    else {
      if(!is.null(nma$temp))
        f <- f*kb*nma$temp
    }
    
    return(f)
  }
