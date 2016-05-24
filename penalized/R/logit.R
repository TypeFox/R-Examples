.logitfit <- function(response, offset) {

  # Finds local gradient and subject weights
  fit <- function(lp, leftout) {
                       
    if (!missing(leftout)) {
      response <- response[!leftout]
      offset <- offset[!leftout]
    }
    
    lp0 <- lp
    if (!is.null(offset)) lp <- lp + offset
    explp <- exp(lp)
    probs <- explp / (1+explp)
    ws <- probs * (1-probs)

    # The residuals
    residuals <- response - probs

    # The loglikelihood
    # loglik <- sum(log(probs[response == 1])) + sum(log(1-probs[response == 0]))
    loglik <- sum(lp[response==1] - log(1+explp[response==1])) + sum(-lp[response==0] - log(1+1/explp[response==0]))
    
    if (!is.na(loglik) && (loglik == -Inf)) loglik <- NA

    return(list(residuals = residuals, loglik = loglik, W = ws, lp = lp, lp0 = lp0, fitted = probs, nuisance = list()))
  }

  # cross-validated likelihood
  cvl <- function(lp, leftout) {
    if (!is.null(offset)) lp <- lp + offset
    probs <- exp(lp) / (1+exp(lp))
    return(sum(log(probs[response == 1 & leftout])) + sum(log(1-probs[response == 0 & leftout])))
  }

  # cross-validated prediction
  prediction <- function(lp, nuisance, which) {
    if (!is.null(offset)) lp <- lp + offset[which]
    out <- exp(lp) / (1+exp(lp))
    out
  }

  return(list(fit = fit, cvl = cvl, prediction = prediction))
}


# mapping from the linear predictor lp to an actual prediction
.logitpredict <- function(lp, nuisance) {
  out <- exp(lp) / (1+exp(lp))
  out
}


# merges predicted probalities
.logitmerge <- function(predictions, groups) {
  out <- unlist(predictions)[sort.list(sort.list(groups))]
  out
}
