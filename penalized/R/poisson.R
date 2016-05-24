.poissonfit <- function(response, offset) {
 
  # Finds local gradient and subject weights
  fit <- function(lp, leftout) {

    if (!missing(leftout)) {
      response <- response[!leftout]
      offset <- offset[!leftout]
    }
    
    lp0 <- lp
    if (!is.null(offset)) lp <- lp + offset
    lambda <- exp(lp)
    ws <- lambda

    # The residuals
    residuals <- response - lambda

    # The loglikelihood
    loglik <- sum(response * log(lambda)) - sum(lambda) - sum(lfactorial(response))
    if (!is.na(loglik) && (loglik == - Inf)) loglik <- NA

    return(list(residuals = residuals, loglik = loglik, W = ws, lp = lp, lp0 = lp0, fitted = lambda, nuisance = list()))
  }

  cvl <- function(lp, leftout) {
    if (!is.null(offset)) lp <- lp + offset
    lambda <- exp(lp[leftout])
    respl <- response[leftout]
    return(sum(respl * log(lambda)) - sum(lambda) - sum(lfactorial(respl)))
  }

  # crossvalidated prediction
  prediction <- function(lp, nuisance, which) {
    if (!is.null(offset)) lp <- lp + offset[which]
    out <- exp(lp)
    out
  }

  return(list(fit = fit, cvl = cvl, prediction = prediction))
}


# mapping from the linear predictor lp to an actual prediction
.poissonpredict <- function(lp, nuisance) {
  out <- exp(lp)
  out
}


# merges predicted probalities
.poissonmerge <- function(predictions, groups) {
  out <- unlist(predictions)[sort.list(sort.list(groups))]
  out
}
