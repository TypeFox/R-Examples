iterate <- function(object, moresteps=200, minlam=0, verbose=FALSE) {
  cl = match.call()

  # Cannot iterate if path completed
  if (object$completepath) stop("Path has completed, nothing to iterate!")

  # Compute the maximum number of steps
  maxsteps = length(object$lambda) + moresteps

  if (object$pathobjs$type %in% c("tall","wide","wide.sparse")) {
    # Iterate dualpathTall
    if (object$pathobjs$type == "tall") {
      out = dualpathTall(object=object,maxsteps=maxsteps,minlam=minlam,verbose=verbose)
    }
    # Iterate dualpathWide
    if (object$pathobjs$type == "wide") {
      out = dualpathWide(object=object,maxsteps=maxsteps,minlam=minlam,verbose=verbose)
    }
    # Iterate dualpathWideSparse
    if (object$pathobjs$type == "wide.sparse") {
      out = dualpathWideSparse(object=object,maxsteps=maxsteps,minlam=minlam,verbose=verbose)
    }

    # Create output object (work usually done by dualpath or genlasso)
    lambda = n0 = y0 = j = D0 = coldif = Xi = X = y2 = NULL
    for (i in 1:length(object)) {
      if (names(object)[i] != "pathobjs") {
        assign(names(object)[i], object[[i]])
      }
    }
    for (i in 1:length(object$pathobjs)) {
      assign(names(object$pathobjs)[i], object$pathobjs[[i]])
    }
    lams = lambda

    # Now we (potentially) need to fix beta, fit, y, bls
    # (This would have been done by dualpath)
    out$df = out$df + coldif
    beta = matrix(y0,n0,length(out$lambda))
    beta[j,] = as.matrix(y0[j] - t(D0[,j])%*%out$u)
    colnames(beta) = colnames(out$u)
    out$beta = beta
    out$fit = beta
    out$y = y0
    out$bls = y0

    out$pathobjs$n0 = n0
    out$pathobjs$y0 = y0
    out$pathobjs$j  = j
    out$pathobjs$D0 = D0
    out$pathobjs$coldif = coldif

    # Now we (potentially) need to account for an X matrix
    # (This would have been done by genlasso)
    if (!is.null(object$X)) {
      out$pathobjs$y2 = y2
      out$pathobjs$Xi = Xi

      out$beta = Xi %*% out$fit
      out$fit = X %*% out$beta
      out$y = object$y
      out$bls = Xi %*% y2
      out$X = X
    }

    # Finally, check if we need to put in the trend order
    # and underlying positions
    if (!is.null(object$ord)) out$ord = object$ord
    if (!is.null(object$pos)) out$pos = object$pos

    out$call = cl
    class(out) = class(object)
    return(out)
  }

  else if (object$pathobjs$type %in% c("fused","fused.l1","fused.x","fused.l1.x")) {
    # Iterate for dualpathFused
    if (object$pathobjs$type == "fused") {
      out = dualpathFused(object=object,maxsteps=maxsteps,minlam=minlam,verbose=verbose)
    }
    # Iterate for dualpathFusedL1
    if (object$pathobjs$type == "fused.l1") {
      out = dualpathFusedL1(object=object,maxsteps=maxsteps,minlam=minlam,verbose=verbose)
    }
    # Iterate for dualpathFusedX
    if (object$pathobjs$type == "fused.x") {
      out = dualpathFusedX(object=object,maxsteps=maxsteps,minlam=minlam,verbose=verbose)
    }
    # Iterate for dualpathFusedL1X
    if (object$pathobjs$type == "fused.l1.x") {
      out = dualpathFusedL1X(object=object,maxsteps=maxsteps,minlam=minlam,verbose=verbose)
    }

    # May need to put in underlying positions (only for
    # the 1d fused lasso)
    if (!is.null(object$pos)) out$pos = object$pos

    out$call = cl
    class(out) = class(object)
    return(out)
  }

  else if (object$pathobjs$type == "trend.x") {
    out = dualpathTrendX(object=object,maxsteps=maxsteps,minlam=minlam,verbose=verbose)
    out$ord = object$ord

    out$call = cl
    class(out) = class(object)
    return(out)
  }

  else {
    # Not a dualpath type that we recognize
    stop(paste("Cannot iterate on path of type",object$pathobjs$type))
  }
}
