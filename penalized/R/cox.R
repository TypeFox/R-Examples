.coxfit <- function(response, offset, strata) {

  n <- nrow(response)
  type <- attr(response, "type")
  if (!type %in% c("right", "counting"))
    stop("Cox model doesn't support \"", type, "\" survival data")

  if (ncol(response) == 2) {
    time <- response[,1]
    status <- response[,2]
    dtimes <- time[status == 1]
    Riskset <- outer(time, dtimes, ">=")
  } else {
    time <- response[,2]
    start <- response[,1]
    status <- response[,3]
    dtimes <- time[status==1]
    Riskset <- outer(time, dtimes, ">=") & outer(start, dtimes, "<")
  }
  whichd <- which(status==1)
  if (!is.null(strata)) {
    dstrata <- strata[status==1]
    Riskset <- Riskset & outer(strata, dstrata, "==")
  } else {
    strata <- factor(rep("baseline", nrow(response)))
    dstrata <- strata[status==1]
  }
  

  # Finds local gradient and subject weights
  fit <- function(lp, leftout) {

    if (!missing(leftout)) {
      status <- status[!leftout]
      time <- time[!leftout]
      strata <- strata[!leftout]
      dleftout <- leftout[whichd]
      dtimes <- dtimes[!dleftout]
      dstrata <- dstrata[!dleftout]
      Riskset <- Riskset[!leftout, !dleftout]
      offset <- offset[!leftout]
    }
    
    lp0 <- lp
    if (!is.null(offset)) lp <- lp + offset
    ws <- drop(exp(lp))
    if (any(ws == Inf | ws == 0)) {
      ws <- 1e-10 + 1e10 * status
      exploded <- TRUE
    } else {
      exploded <- FALSE
    }

    breslows <- drop(1 / ws %*% Riskset)
    breslow <- drop(Riskset %*% breslows)

    # The martingale residuals
    residuals <- status - breslow * ws

    # The loglikelihood
    if (!exploded)
      loglik <- -sum(ws * breslow) + sum(log(breslows)) + sum(lp[status==1])
    else
      loglik <- NA
    if (is.na(loglik) || loglik==-Inf) loglik <- NA

    # The weights matrix
    Pij <- outer(ws, breslows) 
    Pij[!Riskset] <- 0
    W <- list(P = Pij, diagW = breslow * ws)        # construct: W = diag(diagW) - P %*% t(P)

    # The fitted baseline(s)
    baseline <- function() {
      sortlistdtimes <- sort.list(dtimes)
      baselines <- lapply(levels(strata), function(stratum) {
        stratumdtimes <- sortlistdtimes[dstrata == stratum]
        if (length(stratumdtimes) > 0) {
          sdtimes <- dtimes[stratumdtimes]
          basecumhaz <- cumsum(breslows[stratumdtimes])
          uniquetimes <- c(TRUE, as.logical(sapply(seq_along(sdtimes)[-length(sdtimes)], function(i) sdtimes[i] != sdtimes[i+1])))
          basesurv <- exp(-basecumhaz[uniquetimes])
          if (max(dtimes) < max(time)) {
            basetimes <- c(sdtimes[uniquetimes], max(time))
            basesurv <- c(basesurv, basesurv[length(basesurv)])
          } else {
            basetimes <- c(sdtimes[uniquetimes])
            basesurv <- c(basesurv)
          }
          if (min(time) > 0) {
            basetimes <- c(0, basetimes)
            basesurv <- c(1, basesurv)
          }
        } else {
          basetimes <- 0
          basesurv <- 1
        }
    
        baseline <- new("breslow")
        baseline@time <- basetimes
        baseline@curves <- matrix(basesurv,1,byrow=TRUE)
        baseline
      })
      groups <- seq_along(levels(strata))
      names(groups) <- levels(strata)
      .coxmerge(baselines, groups)
    }

    return(list(residuals = residuals, loglik = loglik, W = W, lp = lp, lp0 = lp0, fitted = exp(lp), nuisance = list(baseline = baseline)))
  }

  #cross-validated likelihood
  cvl <- function(lp, leftout)
  {
    if (!is.null(offset)) lp <- lp + offset
    ws <- exp(lp)
    somw <- apply(Riskset, 2, function(rr) sum(ws[rr]))
    cvls <- numeric(length(leftout))
    if (sum(leftout) == 1) { # leave-one-out
      for (k in which(leftout)) {
        pij <- ws[k] / somw
        if (status[k] == 1) {
          dk <- which(whichd == k)
          cvls[k] <- sum(log(1 - pij[Riskset[k,] & (seq_along(dtimes) != dk)])) + log(pij[dk])
        } else {
          cvls[k] <- sum(log(1 - pij[Riskset[k,]]))
        }
      }
      return(sum(cvls[leftout]))
    } else {  # k-fold
      PLall <- sum(log(ws[whichd] / somw))
      newsomw <- apply(Riskset, 2, function(rr) sum(ws[rr & !leftout]))
      newwhichd <- whichd[!(whichd %in% which(leftout))]
      PLrest <- sum(log(ws[newwhichd] / newsomw[!(whichd %in% which(leftout))]))
      return(PLall-PLrest)
    }
  }
  
  # cross-validated predictions
  prediction <- function(lp, nuisance, which) {
    if (!is.null(offset)) lp <- lp + offset[which]
    out <- nuisance$baseline()
    out@curves <- out@curves[strata[which],,drop=FALSE]
    out@curves <- out@curves ^ matrix(exp(lp), nrow(out@curves), ncol(out@curves))
    out
  }

  return(list(fit = fit, cvl = cvl, prediction = prediction))
}


# mapping from the linear predictor lp to an actual prediction
.coxpredict <- function(lp, nuisance, strata) {
  if (is.null(strata)) strata <- rep(1, length(lp))
  out <- nuisance$baseline
  out@curves <- out@curves[strata,,drop=FALSE]
  out@curves <- out@curves ^ matrix(exp(lp), nrow(out@curves), ncol(out@curves))
  row.names(out@curves) <- names(lp)
  out
}


# merges predicted survival curves with different time points
# input: a list of breslow objects
.coxmerge <- function(predictions, groups) {
  times <- sort(unique(unlist(lapply(predictions, time))))
  curves <- lapply(predictions, function(pred) {
    res <- matrix(NA, nrow(pred@curves), length(times))
    res[,times %in% time(pred)] <- pred@curves
    # We interpolate all NAs except in the tail
    startnas <- is.na(res[,1])
    if (any(startnas)) res[startnas,1] <- 1
    endnas <- rev(cumsum(is.na(rev(res[1,])))==1:ncol(res))
    ready <- !any(is.na(res[1,!endnas]))
    while (!ready) {
      nas <- which(is.na(res[1,!endnas]))
      ready <- !any(is.na(res[1,nas-1]))
      res[,nas] <- res[,nas-1]
    }
    res
  })
  out <- new("breslow")
  out@time <- times
  out@curves <- matrix(0, sum(sapply(curves, nrow)), length(times))
  for (i in 1:length(curves)) {
    out@curves[groups==i,] <- curves[[i]]
  }
  rownames(out@curves) <- names(groups)
  out
}


