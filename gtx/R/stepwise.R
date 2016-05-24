# specification of models assumes that real biggest model is union of
# biggest and smallest arguments, i.e. no need to specify terms twice

stepup.moments2 <- function(xtwx, leftvar, biggest, smallest, p.thresh = 0.05, n = NULL, vscale = NULL) {
  ## estimating vscale from data assumes normal/identity model
  ## if using vscale!=1 check definition vis-a-vis standard GLM notation, maybe inverse?
  stopifnot(is.matrix(xtwx))
  stopifnot(nrow(xtwx) == ncol(xtwx))
  stopifnot(all(rownames(xtwx) == colnames(xtwx)))
  if (is.null(n)) stopifnot("ONE" %in% rownames(xtwx))
  stopifnot(length(leftvar) == 1)
  stopifnot(leftvar %in% rownames(xtwx))
  stopifnot(length(biggest) >= 1)
  stopifnot(all(biggest %in% rownames(xtwx)))
  stopifnot(length(smallest) >= 1)
  stopifnot(all(smallest %in% rownames(xtwx)))
  n.arg <- n # cache this for calling est.moments2 with original calling arguments at end
  vscale.arg <- vscale # ditto
  
  if (is.null(n)) n <- xtwx["ONE", "ONE"]

  lidx <- match(leftvar, rownames(xtwx)) # LHS of model
  sidx <- match(smallest, rownames(xtwx)) # smallest model
  bidx <- match(biggest, rownames(xtwx)) # biggest model
  
  nidx <- sidx # initialise the current "null" model
  myxtwxi <- try(solve(xtwx[nidx, nidx, drop = FALSE]), silent = TRUE)
  if (class(myxtwxi) == "matrix") { # X'X is nonsingular
    myxty <- xtwx[nidx, lidx, drop = FALSE]
    if (is.null(vscale)) {
      ssr <- as.double(xtwx[lidx, lidx, drop = TRUE] - t(myxty) %*% myxtwxi %*% myxty)
      vscale <- ssr/(n - length(nidx))
    }
    ll.nidx <- as.double(t(myxty) %*% myxtwxi %*% myxty) / vscale
  } else {
    stop("Smallest model has non-identifiable coefficients")
  }

  done <- FALSE
  chi2.thresh <- qchisq(p.thresh, df = 1, lower.tail = FALSE)
  
  while(!done) {
    delta <- rep(NA, nrow(xtwx)) # store log-likelihood improvement in model fit
    for (tidx in setdiff(bidx, nidx)) { # tidx is "testing" variable
      widx <- c(nidx, tidx) # working model includes tidx
      myxtwxi <- try(solve(xtwx[widx, widx, drop = FALSE]), silent = TRUE)
      if (class(myxtwxi) == "matrix") { # X'X is nonsingular
        myxty <- xtwx[widx, lidx, drop = FALSE]
        if (is.null(vscale)) {
          ssr <- as.double(xtwx[lidx, lidx, drop = TRUE] - t(myxty) %*% myxtwxi %*% myxty)
          vscale <- ssr/(n - length(widx))
        }
        delta[tidx] <- as.double(t(myxty) %*% myxtwxi %*% myxty) / vscale - ll.nidx
      } # else X'X is singular, assume no improvement in model fit
    }
    if (any(!is.na(delta)) && max(delta, na.rm = TRUE) > chi2.thresh) {
      aidx <- which.max(delta)
      cat("Adding", rownames(xtwx)[aidx], "P-value =", signif(pchisq(delta[aidx], df = 1, lower.tail = FALSE), 3), "\n")
      nidx <- c(nidx, aidx)
      ll.nidx <- ll.nidx + delta[aidx]
    } else {
      done <- TRUE
    }
  }
  
  return(est.moments2(xtwx, leftvar, rownames(xtwx)[nidx], n = n.arg, vscale = vscale.arg))
}

stepdown.moments2 <- function(xtwx, leftvar, biggest, smallest, p.thresh = 0.05, n = NULL, vscale = NULL) {
  ## estimating vscale from data assumes normal/identity model
  ## if using vscale!=1 check definition vis-a-vis standard GLM notation, maybe inverse?
  stopifnot(is.matrix(xtwx))
  stopifnot(nrow(xtwx) == ncol(xtwx))
  stopifnot(all(rownames(xtwx) == colnames(xtwx)))
  if (is.null(n)) stopifnot("ONE" %in% rownames(xtwx))
  stopifnot(length(leftvar) == 1)
  stopifnot(leftvar %in% rownames(xtwx))
  stopifnot(length(biggest) >= 1)
  stopifnot(all(biggest %in% rownames(xtwx)))
  stopifnot(length(smallest) >= 1)
  stopifnot(all(smallest %in% rownames(xtwx)))
  n.arg <- n # cache this for calling est.moments2 with original calling arguments at end
  vscale.arg <- vscale # ditto
  
  if (is.null(n)) n <- xtwx["ONE", "ONE"]

  lidx <- match(leftvar, rownames(xtwx)) # LHS of model
  sidx <- match(smallest, rownames(xtwx)) # smallest model
  bidx <- match(biggest, rownames(xtwx)) # biggest model
  
  nidx <- union(sidx, bidx) # initialise at the actual biggest model
  myxtwxi <- try(solve(xtwx[nidx, nidx, drop = FALSE]), silent = TRUE)
  if (class(myxtwxi) == "matrix") { # X'X is nonsingular
    myxty <- xtwx[nidx, lidx, drop = FALSE]
    if (is.null(vscale)) {
      ssr <- as.double(xtwx[lidx, lidx, drop = TRUE] - t(myxty) %*% myxtwxi %*% myxty)
      vscale <- ssr/(n - length(nidx))
    }
    ll.nidx <- as.double(t(myxty) %*% myxtwxi %*% myxty) / vscale
  } else {
    stop("Biggest model has non-identifiable coefficients")
  }

  done <- FALSE
  chi2.thresh <- qchisq(p.thresh, df = 1, lower.tail = FALSE)
  
  while(!done) {
    delta <- rep(NA, nrow(xtwx)) # store log-likelihood improvement in model fit
    for (tidx in setdiff(nidx, sidx)) { # tidx is "testing" variable
      widx <- setdiff(nidx, tidx) # working model excludes tidx
      myxtwxi <- try(solve(xtwx[widx, widx, drop = FALSE]), silent = TRUE)
      if (class(myxtwxi) == "matrix") { # X'X is nonsingular
        myxty <- xtwx[widx, lidx, drop = FALSE]
        if (is.null(vscale)) {
          ssr <- as.double(xtwx[lidx, lidx, drop = TRUE] - t(myxty) %*% myxtwxi %*% myxty)
          vscale <- ssr/(n - length(widx))
        }
        delta[tidx] <- ll.nidx - as.double(t(myxty) %*% myxtwxi %*% myxty) / vscale
      } # else X'X is singular, this should not happen!!!
    }
    if (any(!is.na(delta)) && min(delta, na.rm = TRUE) <= chi2.thresh) {
      aidx <- which.min(delta)
      cat("Dropping", rownames(xtwx)[aidx], "P-value =", signif(pchisq(delta[aidx], df = 1, lower.tail = FALSE), 3), "\n")
      nidx <- setdiff(nidx, aidx)
      ll.nidx <- ll.nidx - delta[aidx]
    } else {
      done <- TRUE
    }
  }
  
  return(est.moments2(xtwx, leftvar, rownames(xtwx)[nidx], n = n.arg, vscale = vscale.arg))
}
