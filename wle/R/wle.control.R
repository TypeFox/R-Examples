#############################################################
#                                                           #
#	wle.glm.control function                            #
#       wle.lm.control function                             #
#       wle.sur.control function                            #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: March, 12, 2010                               #
#	Version: 0.3                                        #
#                                                           #
#	Copyright (C) 2010 Claudio Agostinelli              #
#                                                           #
#############################################################

wle.glm.control <- function(boot=30, group=NULL, num.sol=1, raf=c('GKL', 'PWD', 'HD', 'NED', 'SCHI2'), tau=0.1, cutpoint=0, powerdown=1, delta=NULL, smooth=NULL, asy.smooth=0.031, tol=10^(-6), equal=10^(-3), max.iter=500, window.size=NULL, use.asymptotic=NULL, use.smooth=TRUE, mle.dispersion=FALSE, verbose=FALSE) {  
  if (!is.numeric(boot) || boot <= 0) 
    stop("value of 'boot' must be > 0")
  if(!is.null(group))
    if(!is.numeric(group) || group <= 0)
      stop("value of 'group' must be > 0 or NULL")
  if (!is.numeric(num.sol) || num.sol < 1 || !isTRUE(all.equal(as.integer(num.sol), num.sol)))
    stop("value of 'num.sol' must be an integer")
  raf <- match.arg(raf)
  if ((!is.numeric(tau) || tau > 1 || tau < 0) && raf=='GKL')
    stop("value of 'tau' must be in [0,1] when 'raf=GKL'")
  ##park+basu+2003.pdf
  if ((!is.numeric(tau) || tau < -1) && raf=='PWD')
    stop("value of 'tau' must be in [-1,Inf] when 'raf=PWD'")
  if (!is.numeric(cutpoint) || cutpoint > 1 || cutpoint < 0)
    stop("value of 'cutpoint' must be in [0,1]")
  if (!is.numeric(powerdown))
    stop("value of 'powerdown' must be numeric")  
  ##lindsay+1994.pdf
  if(!is.null(delta))
    if(!is.numeric(delta) || delta <= 0 || delta >=1)
      stop("value of 'delta' must be > 0 and < 1")
  if(!is.null(smooth))
    if(!is.numeric(smooth) || smooth <= 0)
      stop("value of 'smooth' must be > 0 or NULL")
  if(!is.numeric(asy.smooth) || asy.smooth <= 0)
    stop("value of 'asy.smooth' must be > 0")
  if (!is.numeric(tol) || tol <= 0)
    stop("value of 'tol' must be > 0")
  if (!is.numeric(equal) || equal <= tol)
    stop("value of 'equal' must be greater then 'tol'")
  if (!is.numeric(max.iter) || max.iter <= 0) 
    stop("maximum number of iterations must be > 0")
  if(!is.null(window.size))
    if(!is.numeric(window.size) || window.size <= 0)
      stop("value of 'window.size' must be > 0 or NULL")
  if(!is.null(use.asymptotic))
    if(!is.numeric(use.asymptotic) || length(use.asymptotic)!=1 || use.asymptotic <= 0)
      stop("value of 'use.asymptotic' must be a scalar > 0")
  if (!(is.logical(use.smooth) && length(use.smooth)==1))
    stop("'use.smooth' must be logical")
  if (!(is.logical(mle.dispersion) && length(mle.dispersion)==1))
    stop("'mle.dispersion' must be logical")  
  if (!(is.logical(verbose) && length(verbose)==1))
    stop("'verbose' must be logical")  
  list(boot=boot, group=group, num.sol=num.sol, raf=raf, tau=tau, cutpoint=cutpoint, powerdown=powerdown, delta=delta, smooth=smooth, asy.smooth=asy.smooth, tol=tol, equal=equal, max.iter=max.iter, window.size=window.size, use.asymptotic=use.asymptotic, use.smooth=use.smooth, mle.dispersion=mle.dispersion, verbose=verbose)
}

#############################################################
#                                                           #
#	wle.lm.control function                             #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: May, 18, 2010                                 #
#	Version: 0.1                                        #
#                                                           #
#	Copyright (C) 2010 Claudio Agostinelli              #
#                                                           #
#############################################################

wle.lm.control <- function(nstart=30, group=NULL, num.sol=1, raf=c('HD', 'NED', 'SCHI2', 'GKL', 'PWD'), tau=0.1, cutpoint=0, powerdown=1, smooth=0.031, tol=10^(-6), equal=10^(-3), max.iter=500, verbose=FALSE) {
  
  if (!is.numeric(nstart) || nstart <= 0) 
    stop("value of 'nstart' must be > 0")
  if(!is.null(group))
    if(!is.numeric(group) || group <= 0)
      stop("value of 'group' must be > 0 or NULL")
  if (!is.numeric(num.sol) || num.sol < 1 || !isTRUE(all.equal(as.integer(num.sol), num.sol)))
    stop("value of 'num.sol' must be an integer")
  raf <- match.arg(raf)
  if (raf=='GKL' | raf=='PWD')
    stop("no implemented for GKL or PWD RAF for now!")
  if ((!is.numeric(tau) || tau > 1 || tau < 0) && raf=='GKL')
    stop("value of 'tau' must be in [0,1] when 'raf=GKL'")
  ##park+basu+2003.pdf
  if ((!is.numeric(tau) || tau < -1) && raf=='PWD')
    stop("value of 'tau' must be in [-1,Inf] when 'raf=PWD'")
  if (!is.numeric(cutpoint) || cutpoint > 1 || cutpoint < 0)
    stop("value of 'cutpoint' must be in [0,1]")
  if (!is.numeric(powerdown))
    stop("value of 'powerdown' must be numeric")  
  ##lindsay+1994.pdf
  if(!is.numeric(smooth) || smooth <= 0)
    stop("value of 'smooth' must be > 0")
  if (!is.numeric(tol) || tol <= 0)
    stop("value of 'tol' must be > 0")
  if (!is.numeric(equal) || equal <= tol)
    stop("value of 'equal' must be greater then 'tol'")
  if (!is.numeric(max.iter) || max.iter <= 0) 
    stop("maximum number of iterations must be > 0")
  if (!(is.logical(verbose) && length(verbose)==1))
    stop("'verbose' must be logical")  
  list(nstart=nstart, group=group, num.sol=num.sol, raf=raf, tau=tau, cutpoint=cutpoint, powerdown=powerdown, smooth=smooth, tol=tol, equal=equal, max.iter=max.iter, verbose=verbose)
}

#############################################################
#                                                           #
#	wle.sur.control function                            #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: October, 12, 2011                             #
#	Version: 0.1                                        #
#                                                           #
#	Copyright (C) 2011 Claudio Agostinelli              #
#                                                           #
#############################################################

wle.sur.control <- function(first=list(list(nstart=30, group=NULL, num.sol=1, raf=c('HD', 'NED', 'SCHI2', 'GKL', 'PWD'), tau=0.1, cutpoint=0, powerdown=1, smooth=0.031, tol=10^(-6), equal=10^(-3), max.iter=500, verbose=FALSE)), second=list(nstart=30, group=NULL, num.sol=1, raf=c('HD', 'NED', 'SCHI2', 'GKL', 'PWD'), tau=0.1, cutpoint=0, powerdown=1, smooth=0.031, tol=10^(-6), equal=10^(-3), max.iter=500, verbose=FALSE), num.formula=NULL) {
  default <- list(nstart=30, group=NULL, num.sol=1, raf='HD', tau=0.1, cutpoint=0, powerdown=1, smooth=0.031, tol=10^(-6), equal=10^(-3), max.iter=500, verbose=FALSE)
  ### First step
  if (is.list(first)) {
    nfirst <- length(first)
    for (i in 1:nfirst) {
      dn <- setdiff(names(default), names(first[[i]]))
      first[[i]] <- c(first[[i]], default[dn])
      first[[i]]$raf <- first[[i]]$raf[1]
      first[[i]]$raf <- match.arg(first[[i]]$raf, c('HD', 'NED', 'SCHI2', 'GKL', 'PWD'))
    }
    if (is.null(num.formula)) {
      num.formula <- 2
    } else {
      if (nfirst !=1 && nfirst != num.formula)
        stop("the length of 'first', when different from 1, must be equal to 'num.formula' when specified")
    }
    num.formula <- as.integer(num.formula)
    if (is.na(num.formula) || !is.numeric(num.formula) || num.formula < 2)
      stop("'num.formula' must be an a positive integer not less than 2")
    if (nfirst == 1) {
      for (i in 2:num.formula) {
        first[[i]] <- first[[1]]
      }
    } 
    post <- switch(nfirst, "st", "nd", "rd")
    if (is.null(post))
      post <- "th"
    for (i in 1:num.formula) {
      if (!is.numeric(first[[i]]$nstart) || first[[i]]$nstart <= 0) 
        stop("value of 'nstart' must be > 0 in ", i, post, " set of control parameters")
      if(!is.null(first[[i]]$group))
        if(!is.numeric(first[[i]]$group) || first[[i]]$group <= 0)
          stop("value of 'group' must be > 0 or NULL in ", i, post, " set of control parameters")
      if (!is.numeric(first[[i]]$num.sol) || first[[i]]$num.sol < 1 || !isTRUE(all.equal(as.integer(first[[i]]$num.sol), first[[i]]$num.sol)))
        stop("value of 'num.sol' must be an integer in ", i, post, " set of control parameters")
      if (first[[i]]$raf=='GKL' | first[[i]]$raf=='PWD')
        stop("no implemented for GKL or PWD RAF for now in ", i, post, " set of control parameters")
      if ((!is.numeric(first[[i]]$tau) || first[[i]]$tau > 1 || first[[i]]$tau < 0) && first[[i]]$raf=='GKL')
        stop("value of 'tau' must be in [0,1] when 'raf=GKL' in ", i, post, " set of control parameters")
  ##park+basu+2003.pdf
      if ((!is.numeric(first[[i]]$tau) || first[[i]]$tau < -1) && first[[i]]$raf=='PWD')
        stop("value of 'tau' must be in [-1,Inf] when 'raf=PWD' in ", i, post, " set of control parameters")
      if (!is.numeric(first[[i]]$cutpoint) || first[[i]]$cutpoint > 1 || first[[i]]$cutpoint < 0)
        stop("value of 'cutpoint' must be in [0,1] in ", i, post, " set of control parameters")
      if (!is.numeric(first[[i]]$powerdown))
        stop("value of 'powerdown' must be numeric in ", i, post, " set of control parameters")
  ##lindsay+1994.pdf
      if(!is.numeric(first[[i]]$smooth) || first[[i]]$smooth <= 0)
        stop("value of 'smooth' must be > 0 in ", i, post, " set of control parameters")
      if (!is.numeric(first[[i]]$tol) || first[[i]]$tol <= 0)
        stop("value of 'tol' must be > 0 in ", i, post, " set of control parameters")
      if (!is.numeric(first[[i]]$equal) || first[[i]]$equal <= first[[i]]$tol)
        stop("value of 'equal' must be greater then 'tol' in ", i, post, " set of control parameters")
      if (!is.numeric(first[[i]]$max.iter) || first[[i]]$max.iter <= 0) 
        stop("maximum number of iterations must be > 0 in ", i, post, " set of control parameters")
      if (!(is.logical(first[[i]]$verbose) && length(first[[i]]$verbose)==1))
        stop("'verbose' must be logical in ", i, post, " set of control parameters")
    }
  }
### Second step
  dn <- setdiff(names(default), names(second))
  second <- c(second, default[dn])
  if (!is.numeric(second$nstart) || second$nstart <= 0) 
    stop("value of 'nstart' must be > 0 in the second step")
  if (!is.null(second$group))
    if (!is.numeric(second$group) || second$group <= 0)
      stop("value of 'group' must be > 0 or NULL in the second step")
    if (!is.numeric(second$num.sol) || second$num.sol < 1 || !isTRUE(all.equal(as.integer(second$num.sol), second$num.sol)))
      stop("value of 'num.sol' must be an integer in the second step")
    second$raf <- second$raf[1]
    second$raf <- match.arg(second$raf, c('HD', 'NED', 'SCHI2', 'GKL', 'PWD'))
    if (second$raf=='GKL' | second$raf=='PWD')
      stop("no implemented for GKL or PWD RAF for now in the second step")
    if ((!is.numeric(second$tau) || second$tau > 1 || second$tau < 0) && second$raf=='GKL')
      stop("value of 'tau' must be in [0,1] when 'raf=GKL' in the second step")
  ##park+basu+2003.pdf
    if ((!is.numeric(second$tau) || second$tau < -1) && second$raf=='PWD')
      stop("value of 'tau' must be in [-1,Inf] when 'raf=PWD' in the second step")
    if (!is.numeric(second$cutpoint) || second$cutpoint > 1 || second$cutpoint < 0)
      stop("value of 'cutpoint' must be in [0,1] in the second step")
    if (!is.numeric(second$powerdown))
      stop("value of 'powerdown' must be numeric in the second step")
  ##lindsay+1994.pdf
    if (!is.numeric(second$smooth) || second$smooth <= 0)
      stop("value of 'smooth' must be > 0 in the second step")
    if (!is.numeric(second$tol) || second$tol <= 0)
      stop("value of 'tol' must be > 0 in the second step")
    if (!is.numeric(second$equal) || second$equal <= second$tol)
      stop("value of 'equal' must be greater then 'tol' in the second step")
    if (!is.numeric(second$max.iter) || second$max.iter <= 0) 
      stop("maximum number of iterations must be > 0 in the second step")
    if (!(is.logical(second$verbose) && length(second$verbose)==1))
      stop("'verbose' must be logical in the second step")
  result <- list(first=first, second=second)
  return(result)
}




