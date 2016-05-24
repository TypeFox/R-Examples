mlen <-
function(data, coords, loc.id, direction, mle = "avg") {
  # Empirical estimation mean-lengths (for embeded data)
  #
  #      data vector of data or 
  #    coords matrix of coordinates
  #    loc.id location Id (which_lines output)
  # direction vector (orres versor) of choosen direction
  #       mle string of characters, 
  #           if "mlk" the MLEs will be returned (log-normal distro assumed)
  #           it can be logical for backward compatibility reasons

  mle <- mle[1L]
  if (is.logical(mle)) {
    mle <- ifelse(mle, "mlk", "avg") 
  }
  else {
    if (is.character(mle)) {
      if (!mle %in% c("avg", "mlk", "trm", "mdn")) mle <- "avg"      
    }
    else {
      mle <- "avg"
    }
  }
  # Mean-Length Estimation via method of moments (averaging)
  if (mle == "avg") {
    gl <- getlen(data, coords, loc.id, direction, zero.allowed = TRUE)
    meanlen <- tapply(gl$length + gl$maxcens, gl$categories, mean)
  }
  # Mean-Length Estimation via maximum likelihood (log-normal distribution)
  if (mle == "mlk") {
    gl <- getlen(data, coords, loc.id, direction, zero.allowed = TRUE)
    nk <- nlevels(data)
    if (length(data) < nk) stop("there are not enough data to estimate the parameters")
    param <- vector("numeric", 2 * nk)
    gl$categories <- as.integer(gl$categories)

    NegLik <- function(param) {
      pus <- plnorm(gl$length + gl$maxcens, meanlog = param[gl$categories],
                    sdlog = exp(param[gl$categories + nk]))
      pls <- plnorm(gl$length, meanlog = param[gl$categories],
                    sdlog = exp(param[gl$categories + nk]))
      return(- sum(log(abs(pus - pls) + .Machine$double.neg.eps)))
    }

    res <- nlminb(param, NegLik, lower = -Inf, upper = Inf)
    message("Optimization message: ", res$message, sep = "")
    meanlen <- exp(res$par[1:nk] + 0.5 * exp(2 * res$par[(nk + 1):(2 * nk)]))
  }
# Mean-Length Estimation via trimmed averaging
  if (mle == "trm") {
    gl <- getlen(data, coords, loc.id, direction, zero.allowed = FALSE)
    meanlen <- tapply(gl$length, gl$categories, mean)
  }
# Mean-Length Estimation via trimmed median calculation
  if (mle == "mdn") {
    gl <- getlen(data, coords, loc.id, direction, zero.allowed = TRUE)
    meanlen <- tapply(gl$length + gl$maxcens, gl$categories, median)
  }
  iiff <- is.finite(meanlen)
  zzzz <- meanlen[iiff] <= 0
  if (!all(iiff) || any(zzzz)) {
    stdlen <- apply(coords, 2, function(x) mean(diff(sort(unique(x)))))
    stdlen <- sqrt(sum(abs(direction) / sqrt(sum(direction^2)) * stdlen^2))
    meanlen[!iiff] <- stdlen
    meanlen[iiff][zzzz] <- stdlen
  }
  return(meanlen)
}
