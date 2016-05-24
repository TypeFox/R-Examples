#  $Id: csurv.R,v 1.6 2003/03/28 12:55:32 hothorn Exp $

csurv <- function(newdata, pred, minprob=0, window=0.0001) {

  N <- nrow(newdata)
  if (!"hazard" %in% names(attributes(newdata)))
    stop("hazards attribute to newdata missing")
  hazards <- attr(newdata, "hazard")

  error <- rep(0, N)

  # if there is only one prediction for all observations
  GETPROB <- TRUE
  if (inherits(pred, "survfit")) {
    times <- pred$time 			# get times
    predprob <- getsurv(pred, times)	# get steps
    GETPROB <- FALSE
  }

  for (i in 1:N) {
    if (GETPROB) {
      times <- pred[[i]]$time 		# get times
      predprob <- getsurv(pred[[i]], times)	# get steps
    }
    # compute the integrated squared difference between
    # KM and S(t)
    # minprob: stop integration when S(t) < minprob
    lasttime <- -(log(minprob) / hazards[i])
    if (max(times) > lasttime) {
      thisprob  <- predprob[times <= lasttime]
      thistimes <- times[times <= lasttime]
    } else {
      thisprob  <- predprob
      thistimes <- times
    }
    error[i] <- .Call("SdiffKM", as.double(c(0,thistimes)), 
                       as.double(c(1,thisprob)),
                       as.double(c(hazards[i], window)), PACKAGE="ipred")
    # adjust for time scale by last event
    error[i] <- error[i]/max(thistimes)
    if (length(unique(hazards)) == 1) {
      error <- error[i]
      break
    }
  }
  error <- mean(error)
  error
}

foo <- function (time, prob, hazard, window) 
{
    myint <- 0
    time <- c(0, time)
    s <- exp(-time * hazard)
    prob <- c(1, prob)
    for (i in 1:(length(time)-1)) {   
        d <- time[i+1] - time[i]
        if (d < window) {
            myint <- myint + 0.5 * d * ((prob[i] - s[i])^2 +
                (prob[i] - s[i + 1])^2)
        }
        else {
            k <- ceiling(d/window)
            wi <- d/k
            for (j in 1:k) myint <- myint + 0.5 * wi * ((prob[i] -
                exp(-(time[i] + (j - 1) * wi) * hazard))^2 +
                (prob[i] - exp(-(time[i] + j * wi) * hazard))^2)
        }
    }
    myint
}
