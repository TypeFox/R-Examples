.lmfit <- function(response, offset) {

  # Finds local gradient and subject weights
  fit <- function(lp, leftout) {
    if (!missing(leftout)) {
      response <- response[!leftout]
      offset <- offset[!leftout]
    }

    # The residuals
    lp0 <- lp
    if (!is.null(offset)) lp <- lp + offset
    residuals <- drop(response - lp)

    # The loglikelihood
    ss <- sum(residuals * residuals)
    if (missing(leftout)) n <- length(lp) else n <- sum(!leftout)
    loglik <- (-n/2) * (log(2*pi/n) + 1 + log(ss + .Machine$double.xmin))

    return(list(residuals = residuals, loglik = loglik, W = 1, lp = lp, lp0 = lp0, fitted = lp, nuisance = list(sigma2 = ss/n)))
  }

  # cross-validated likelihood
  cvl <- function(lp, leftout) {
    if (!is.null(offset)) lp <- lp + offset
    residuals <- response - lp
    sigma2 <- sum(residuals[!leftout] * residuals[!leftout]) / sum(!leftout)
    ss <- sum(residuals[leftout] * residuals[leftout])

    return(-(sum(leftout)/2) * log(2*pi*sigma2) - ss / (2*sigma2))
  }
  
  prediction <- function(lp, nuisance, which) {
    if (!is.null(offset)) lp <- lp + offset[which]
    out <- cbind(mu = lp, sigma2 = nuisance$sigma2)
    out
  }

  return(list(fit = fit, cvl = cvl, prediction = prediction))
}


# mapping from the linear predictor lp to an actual prediction
.lmpredict <- function(lp, nuisance) {
  out <- drop(cbind(mu = lp, sigma2 = nuisance$sigma2))
  out
}
  

# merges predicted means and variances
.lmmerge <- function(predictions, groups) {

  out <- matrix(0, sum(sapply(predictions, nrow)), 2)
  for (i in 1:length(predictions)) {
    out[groups==i,] <- predictions[[i]]
  }
  colnames(out) <- c("mu", "sigma2")
  rownames(out) <- names(groups)
  out
}
