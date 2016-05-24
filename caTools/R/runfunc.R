#===========================================================================#
# caTools - R library                                                       #
# Copyright (C) 2005 Jarek Tuszynski                                        #
# Distributed under GNU General Public License version 3                    #
#===========================================================================#
#source('C:/Projects/R Packages/caTools/R/runfunc.R')

runmean = function(x, k, alg=c("C", "R", "fast", "exact"),
                   endrule=c("mean", "NA", "trim", "keep", "constant", "func"),
                   align = c("center", "left", "right"))
{
  alg     = match.arg(alg)
  endrule = match.arg(endrule)
  align   = match.arg(align)
  dimx = dim(x) # Capture dimension of input array - to be used for formating y
  x = as.vector(x)
  n = length(x)
  if (k<=1) return (x)
  if (k >n) k = n
  k2 = k%/%2

  if (alg=="exact") {
    y <- .C("runmean_exact", x, y = double(n) , as.integer(n), as.integer(k),
            NAOK=TRUE, PACKAGE="caTools")$y
  } else if (alg=="C") {
    y <- .C("runmean", as.double(x), y = double(n), as.integer(n), as.integer(k),
            NAOK=TRUE, PACKAGE="caTools")$y
  } else if (alg=="fast") {
    y <- .C("runmean_lite", as.double(x), y = double(n), as.integer(n), as.integer(k),
            NAOK=TRUE, PACKAGE="caTools")$y
  } else {     # the similar algorithm implemented in R language
      y = double(n)
    k1 = k-k2-1
    y = c( sum(x[1:k]), diff(x,k) ); # find the first sum and the differences from it
    y = cumsum(y)/k                  # apply precomputed differences
    y = c(rep(0,k1), y, rep(0,k2))   # make y the same length as x
    if (endrule=="mean") endrule="func"
  }
  y = EndRule(x, y, k, dimx, endrule, align, mean, na.rm=TRUE)
  return(y)
}

#==============================================================================

runmin = function(x, k, alg=c("C", "R"),
                  endrule=c("min", "NA", "trim", "keep", "constant", "func"),
                  align = c("center", "left", "right"))
{
  alg = match.arg(alg)
  align   = match.arg(align)
  endrule = match.arg(endrule)
  dimx = dim(x)  # Capture dimension of input array - to be used for formating y
  x = as.vector(x)
  n = length(x)
  if (k<=1) return (x)
  if (k >n) k = n

  if (alg=="C") {
    y <- .C("runmin", as.double(x), y = double(n), as.integer(n), as.integer(k),
            NAOK=TRUE, PACKAGE="caTools")$y
  } else { # the similar algorithm implemented in R language
      y = double(n)
    k2 = k%/%2
    k1 = k-k2-1
    a <- y[k1+1] <- min(x[1:k], na.rm=TRUE)
    if (k!=n) for (i in (2+k1):(n-k2)) {
      if (a==y[i-1]) # point leaving the window was the min, so ...
        y[i] = min(x[(i-k1):(i+k2)], na.rm=TRUE) # recalculate min of the window
      else           # min=y[i-1] is still inside the window
        y[i] = min(y[i-1], x[i+k2 ], na.rm=TRUE) # compare it with the new point
      a = x[i-k1]    # point that will be removed from the window next
      if (!is.finite(a)) a=y[i-1]+1 # this will force the 'else' option
    }
    if (endrule=="min") endrule="func"
  }
  y = EndRule(x, y, k, dimx, endrule, align, min, na.rm=TRUE)
  return(y)
}

#==============================================================================

runmax = function(x, k, alg=c("C", "R"),
                  endrule=c("max", "NA", "trim", "keep", "constant", "func"),
                  align = c("center", "left", "right"))
{
  alg     = match.arg(alg)
  endrule = match.arg(endrule)
  align   = match.arg(align)
  dimx = dim(x) # Capture dimension of input array - to be used for formating y
  x = as.vector(x)
  n = length(x)
  k = as.integer(k)
  if (k<=1) return (x)
  if (k >n) k = n
  y = double(n)

  if (alg=="C") {
    y <- .C("runmax", as.double(x), y = double(n) , as.integer(n), as.integer(k),
            NAOK=TRUE, PACKAGE="caTools")$y
  } else { # the same algorithm implemented in R language
      y = double(n)
    k2 = k%/%2
    k1 = k-k2-1
    a <- y[k1+1] <- max(x[1:k], na.rm=TRUE)
    if (k!=n) for (i in (2+k1):(n-k2)) {
      if (a==y[i-1]) # point leaving the window was the max, so ...
        y[i] = max(x[(i-k1):(i+k2)], na.rm=TRUE) # recalculate max of the window
      else           # max=y[i-1] is still inside the window
        y[i] = max(y[i-1], x[i+k2 ], na.rm=TRUE) # compare it with the new point
      a = x[i-k1]    # point that will be removed from the window next
      if (!is.finite(a)) a=y[i-1]+1 # this will force the 'else' option
    }
    if (endrule=="max") endrule="func"
  }
  y = EndRule(x, y, k, dimx, endrule, align, max, na.rm=TRUE)
  return(y)
}

#==============================================================================

runquantile = function(x, k, probs, type=7,
                endrule=c("quantile", "NA", "trim", "keep", "constant", "func"),
                align = c("center", "left", "right"))
{ ## see http://mathworld.wolfram.com/Quantile.html for very clear definition
  ## of different quantile types
  endrule = match.arg(endrule)
  align   = match.arg(align)
  dimx = dim(x) # Capture dimension of input array - to be used for formating y
  yIsVec = is.null(dimx) # original x was a vector
  x    = as.vector(x)
  n    = length(x)
  np   = length(probs)
  k    = as.integer(k)
  type = as.integer(type)
  if (k<=1) return (rep(x,n,np))
  if (k >n) k = n
  if (is.na(type) || (type < 1 | type > 9))
    warning("'type' outside allowed range [1,9]; changing 'type' to ", type<-7)

  y = double(n*np)
  y <- .C("runquantile", as.double(x), y = y , as.integer(n), as.integer(k),
          as.double(probs), as.integer(np),as.integer(type),
          NAOK=TRUE, PACKAGE="caTools")$y
  dim(y) =  c(n,np)

  for (i in 1:np) {   # for each percentile
    yTmp = EndRule(x, y[,i], k, dimx, endrule, align, quantile, probs=probs[i], type=type, na.rm=TRUE)
    if (i==1) {
      if (is.null(dimx)) dimy = length(yTmp) else dimy = dim(yTmp)
      yy = matrix(0,length(yTmp),np)   # initialize output array
    }
    yy[,i] = as.vector(yTmp)
  }
  if (np>1) dim(yy) = c(dimy,np) else dim(yy) = dimy
  return(yy)
}

#==============================================================================

runmad = function(x, k, center = runmed(x,k), constant = 1.4826,
                  endrule=c("mad", "NA", "trim", "keep", "constant", "func"),
                  align = c("center", "left", "right"))
{
  endrule = match.arg(endrule)
  align   = match.arg(align)
  dimx = dim(x) # Capture dimension of input array - to be used for formating y
  x = as.vector(x)
  n = length(x)
  if (k<3) stop("'k' must be larger than 2")
  if (k>n) k = n
  y <- .C("runmad", as.double(x), as.double(center), y = double(n),
          as.integer(n), as.integer(k), NAOK=TRUE, PACKAGE="caTools")$y
  y = EndRule(x, y, k, dimx, endrule, align, mad, constant=1, na.rm=TRUE)
  return(constant*y)
}

#==============================================================================

runsd = function(x, k, center = runmean(x,k),
                 endrule=c("sd", "NA", "trim", "keep", "constant", "func"),
                 align = c("center", "left", "right"))
{
  endrule = match.arg(endrule)
  align   = match.arg(align)
  dimx = dim(x) # Capture dimension of input array - to be used for formating y
  x = as.vector(x)
  n = length(x)
  if (k<3) stop("'k' must be larger than 2")
  if (k>n) k = n
  y <- .C("runsd", as.double(x), as.double(center), y = double(n),
          as.integer(n), as.integer(k), NAOK=TRUE, PACKAGE="caTools")$y
  y = EndRule(x, y, k, dimx, endrule, align, sd, na.rm=TRUE)
  return(y)
}

#==============================================================================

EndRule = function(x, y, k, dimx,
             endrule=c("NA", "trim", "keep", "constant", "func"),
             align = c("center", "left", "right"), Func, ...)
{
  # Function which postprocess results of running windows functions and cast
  # them in to specified format. On input y is equivalent to
  #   y = runFUNC(as.vector(x), k, endrule="func", align="center")

  # === Step 1: inspects inputs and unify format ===
  align   = match.arg(align)
  k = as.integer(k)
  k2 = k%/%2
  if (k2<1) k2 = 1
  yIsVec = is.null(dimx) # original x was a vector -> returned y will be a vector
  if (yIsVec) dimx=c(length(y),1) # x & y will become 2D arrays
  dim(x) <- dimx
  dim(y) <- dimx
  n = nrow(x)
  m = ncol(x)
  if (k>n) k2 = (n-1)%/%2
  k1 = k-k2-1
  if (align=="center" && k==2) align='right'

  # === Step 2: Apply different endrules ===
  if (endrule=="trim") {
    y = y[(k1+1):(n-k2),] # change y dimensions
  } else if (align=="center") {
    idx1 = 1:k1
    idx2 = (n-k2+1):n
    # endrule calculation in R will be skipped for most common case when endrule
    # is default and array was a vector not a matrix
    if (endrule=="NA") {
      y[idx1,] = NA
      y[idx2,] = NA
    } else if (endrule=="keep") {
      y[idx1,] = x[idx1,]
      y[idx2,] = x[idx2,]
    } else if (endrule=="constant") {
      y[idx1,] = y[k1+1+integer(m),]
      y[idx2,] = y[n-k2+integer(m),]
    } else if (endrule=="func" || !yIsVec) {
      for (j in 1:m) {
        for (i in idx1) y[i,j] = Func(x[1:(i+k2),j], ...)
        for (i in idx2) y[i,j] = Func(x[(i-k1):n,j], ...)
      }
    }
  } else if (align=="left") {
    y[1:(n-k1),] = y[(k1+1):n,]
    idx = (n-k+2):n
    if (endrule=="NA") {
      y[idx,] = NA
    } else if (endrule=="keep") {
      y[idx,] = x[idx,]
    } else if (endrule=="constant") {
      y[idx,] = y[n-k+integer(m)+1,]
    } else {
      for (j in 1:m) for (i in idx) y[i,j] = Func(x[i:n,j], ...)
    }
  } else if (align=="right") {
    y[(k2+1):n,] = y[1:(n-k2),]
    idx = 1:(k-1)
    if (endrule=="NA") {
      y[idx,] = NA
    } else if (endrule=="keep") {
      y[idx,] = x[idx,]
    } else if (endrule=="constant") {
      y[idx,] = y[k+integer(m),]
    } else {
      for (j in 1:m) for (i in idx) y[i,j] = Func(x[1:i,j], ...)
    }
  }

  # === Step 4: final casting and return results ===
  if (yIsVec) y = as.vector(y);
  return(y)
}

