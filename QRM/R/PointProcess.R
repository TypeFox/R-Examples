## Copyright (C) 2013 Marius Hofert, Bernhard Pfaff
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.


##
extremalPP <- function(data, threshold = NA, nextremes = NA, ...){
  if(is.timeSeries(data)){
    values <- as.vector(series(data))
    pos <- as.POSIXct(time(data))
    alltimes <- as.integer(julian(pos, ...))
    starttime <- alltimes[1] - 1
    endtime <- alltimes[length(alltimes)]
  } else {
    values <- as.vector(data)
    n <- length(values)
    alltimes <- (1:n) / n
    starttime <- 0
    endtime <- 1
  }
  if(is.na(nextremes) && is.na(threshold))
    stop("Either a threshold or the number of upper extremes must be supplied.")
  if(!is.na(nextremes)){
    threshold <- findthreshold(values, nextremes)}
    times <- alltimes[values > threshold]
    marks <- values[values > threshold] - threshold
    out <- list(times = times, marks = marks, starttime = starttime,
                endtime = endtime, threshold = threshold)
    class(out) <- c("MPP")
    out
}
##
unmark <- function(PP){
  if (class(PP) != "MPP") stop("Not a marked point process")
  PP$marks <- NULL
  class(PP) <- "PP"
  PP
}
##
plot.MPP <- function(x, ...){
  if (class(x) != "MPP") stop("Not a marked point process")
  starttime <- x$starttime
  endtime <- x$endtime
  times <- x$times
  marks <- x$marks
  plot(times, marks, xlim = c(starttime, endtime),
       ylab = "MPP", type = "h", ...)
  return(invisible(list(times = times, marks = marks)))
}
##
plot.PP <- function(x,...){
  if (class(x) != "PP") stop("Not a point process. Did you call unmark()?")
  starttime <- x$starttime
  endtime <- x$endtime
  times <- x$times
  augLength <- length(times) + 1
  count <- 1:augLength
  stepfunc <- stepfun(times, count, f = 0)
  plot.stepfun(stepfunc, xlim = range(starttime, endtime),
               ylim = range(0, augLength), xlab = "Time",
               ylab = "N events", pch = '',
               main = "Time Pattern of Exceedances",...)
  return(invisible(list(steps = stepfunc)))
}
##
fit.POT <- function(PP, markdens = "GPD", ...){
  starttime <- PP$starttime
  endtime <- PP$endtime
  times <- PP$times
  marks <- PP$marks
  span <- endtime - starttime
  par.ests <- length(times) / span
  names(par.ests) <- c("lambda")
  par.ses <- sqrt(par.ests / span)
  names(par.ses) <- names(par.ests)
  names
  ll.max <- -par.ests * span + length(times) * log(par.ests)
  if (class(PP)=="MPP"){
    mark.model <- switch(markdens, GPD = fit.GPD(marks, 0, ...))
    par.ests <- c(par.ests, mark.model$par.ests)
    names(par.ests) <- c("lambda", names(mark.model$par.ests))
    par.ses <- c(par.ses, mark.model$par.ses)
    names(par.ses) <- names(par.ests)
    ll.max <- ll.max + mark.model$ll.max
  }
  list(par.ests = par.ests, par.ses = par.ses, ll.max = ll.max)
}
##
sePP.negloglik <- function(theta, PP, case){
   theta <- abs(theta)
   times <- PP$times
   marks <- PP$marks
   if (class(PP) != "MPP") marks <- 0.0
   endtime <- PP$endtime
   starttime <- PP$starttime
   mu <- theta[1]
   phi <- theta[2]
   voltheta <- theta[-c(1, 2)]
   evol <- volfunction(times,times,marks,voltheta,case)
   lambda.contrib <- mu + phi * evol
   term1 <- sum(log(lambda.contrib))
   gamma <- theta[3]
   delta <- switch(case, 0.0, theta[4], 0.0, theta[5])
   rho <- switch(case, 0.0, 0.0, theta[4], theta[4])
   if(case <= 2){
     terminsum <- (1 - exp(-gamma * (endtime - times))) / gamma
   } else {
     terminsum <- gamma*(1 - (1 + (endtime - times) / gamma)^(-rho)) / rho
   }
   terminsum <- (1 + delta * marks) * terminsum
   term2 <- mu * (endtime - starttime) + phi * sum(terminsum)
   out <- term2 - term1
   out
}
##
seMPP.negloglik <- function(theta, PP, case, markdens){
   theta <- abs(theta)
   times <- PP$times
   marks <- PP$marks
   endtime <- PP$endtime
   starttime <- PP$starttime
   mu <- theta[1]
   phi <- theta[2]
   voltheta <- theta[-c(1, 2, (length(theta) - 2), (length(theta) - 1), length(theta))]
   evol <- volfunction(times, times, marks, voltheta, case)
   lambda.contrib <- mu + phi * evol
   term1 <- sum(log(lambda.contrib))
   xi <- theta[length(theta) - 2]
   beta <- theta[length(theta) - 1]
   alpha <- theta[length(theta)]
   scale <- beta + alpha * evol
   markdensfunc <- switch(markdens, GPD = dGPD)
   lambda.x <- markdensfunc(marks, xi, scale, log = TRUE)
   term3 <- sum(lambda.x)
   gamma <- theta[3]
   delta <- switch(case, 0.0, theta[4], 0.0, theta[5])
   rho <- switch(case, 0.0, 0.0, theta[4], theta[4])
   if(case <= 2){
     terminsum <- (1 - exp(-gamma * (endtime - times))) / gamma
   } else {
     terminsum <- gamma * (1 - (1 + (endtime - times) / gamma)^(-rho)) / rho
   }
   terminsum <- (1 + delta * marks) * terminsum
   term2 <- mu * (endtime - starttime) + phi * sum(terminsum)
   out <- term2 - term1 - term3
   out
}
##
volfunction <- function(anytimes, times, marks, theta, model){
    SEprocExciteFunc(as.vector(anytimes),
                     as.vector(times),
                     as.vector(marks),
                     as.numeric(theta),
                     as.integer(model))
}
##
plot.sePP <- function(x, ...){
  PP <- x$PP
  starttime <- PP$starttime
  endtime <- PP$endtime
  times <- PP$times
  marks <- PP$marks
  theta <- x$par.ests
  case <- x$case
  anytimes <- starttime:endtime
  voltheta <- theta[-c(1, 2)]
  evol <- volfunction(times, times, marks, voltheta, case)
  intensity <- rep(0, length(anytimes))
  intensity[which(anytimes %in% times, arr.ind = TRUE)]  <- theta[1] + theta[2] * evol
  plot(anytimes, intensity, type = "l", xlim = range(starttime, endtime), xlab = "Time", ylab = "Intensity")
  abline(h = length(times) / (endtime - starttime))
  return(invisible(list(times = anytimes, intensity = intensity)))
}
##
fit.sePP <- function(PP, model = c("Hawkes", "ETAS"), mark.influence = TRUE, std.errs = FALSE, ...){
  model <- match.arg(model)
  starttime <- PP$starttime
  endtime <- PP$endtime
  times <- PP$times
  marks <- PP$marks
  span <- endtime-starttime
  rate <- length(times) / span
  if(!((class(PP) =="MPP") & mark.influence)){
    case <- switch(model, Hawkes = 1, ETAS = 3)
    mark.influence=FALSE
  } else {
    case <- switch(model, Hawkes = 2, ETAS = 4)
  }
  theta <- switch(case, c(rate, 0, 0.1), c(rate, 0, 0.1, 0), c(rate, 0, 0.1, 0.1), c(rate, 0, 0.1, 0.1, 0))
  fit <- nlminb(start=theta, objective = sePP.negloglik, PP = PP, case = case, ...)
  par.ests <- fit$par
  par.ests <- abs(par.ests)
  ll.max <- -sePP.negloglik(par.ests,PP,case)
  nms <- c("mu", "phi")
  addnms <- switch(case, c("gamma"), c("gamma", "delta"), c("gamma", "rho"), c("gamma", "rho", "delta"))
  names(par.ests) <- c(nms, addnms)
  if(std.errs){
    hessian <- hessian(sePP.negloglik, par.ests, PP = PP, case = case)
    varcov <- solve(hessian)
    par.ses <- sqrt(diag(varcov))
    names(par.ses) <- names(par.ests)
  } else {
    par.ses <- rep(NA, length(par.ests))
    names(par.ses) <- names(par.ests)
  }
  tstat <- par.ests / par.ses
  out <- list(PP = PP, par.ests = par.ests, par.ses = par.ses, tstat = tstat,
              model = model, mark.model = FALSE,
              mark.influence = mark.influence, case = case, converged =
              fit$message, ll.max = ll.max, fit = fit)
  class(out) <- "sePP"
  out
}
##
fit.seMPP <- function(PP, markdens = "GPD", model = c("Hawkes", "ETAS"), mark.influence = TRUE, predictable = FALSE, std.errs = FALSE, ...){
  if (class(PP) != "MPP") stop("Not marked point process data")
  model <- match.arg(model)
  marks <- PP$marks
  groundmod <- fit.sePP(PP, model, mark.influence, std.errs)
  par.ests <- groundmod$par.ests
  nms <- names(par.ests)
  par.ses <- groundmod$par.ses
  tstat <- groundmod$tstat
  ll.max <- groundmod$ll.max
  converged <- groundmod$converged
  case <- groundmod$case
  mark.model <- switch(markdens, GPD = fit.GPD(marks, 0))
  par.ests <- c(par.ests, mark.model$par.ests)
  names(par.ests) <- c(nms, names(mark.model$par.ests))
  if(std.errs) par.ses <- c(par.ses, mark.model$par.ses)
  ll.max <- ll.max + mark.model$ll.max
  if(!(mark.model$converged)) converged <- FALSE
  if(predictable){
    nms <- names(par.ests)
    theta <- c(par.ests, 0)
    fit <- nlminb(start = theta, objective = seMPP.negloglik, PP = PP, case = case, markdens = markdens, ...)
    par.ests <- fit$par
    par.ests <- abs(par.ests)
    names(par.ests) <- c(nms, "alpha")
    ll.max <- -seMPP.negloglik(par.ests, PP, case, markdens = markdens)
    converged <- fit$message
    if (std.errs){
       hessian <- hessian(seMPP.negloglik, par.ests, PP = PP, case = case, markdens = markdens)
       varcov <- solve(hessian)
       par.ses <- sqrt(diag(varcov))
    }
  }
  if(!std.errs) par.ses <- rep(NA, length(par.ests))
  names(par.ses) <- names(par.ests)
  tstat <- par.ests / par.ses
  out <- list(PP = PP, par.ests = par.ests, par.ses = par.ses, tstat = tstat,
              model = model, mark.model = TRUE, mark.influence =
              mark.influence, case = case, predictable = predictable, converged
              = converged, ll.max = ll.max)
  class(out) <- "sePP"
  out
}
##
stationary.sePP <- function(sePP){
  mark.model <- sePP$mark.model
  mark.influence <- sePP$mark.influence
  if(mark.model){
    if(sePP$predictable) stop("Only implemented for unpredictable marked point processes")
  }
  not.mark.model <- !(mark.model)
  if (not.mark.model & mark.influence) stop("Need auxiliary model for marks")
  par.ests <- sePP$par.ests
  phi <- par.ests[2]
  gamma <- par.ests[3]
  case <- sePP$case
  if(sePP$model == "Hawkes"){
    eta <- phi / gamma
  } else {
    rho <- switch(case, 0.0, 0.0, par.ests[4], par.ests[4])
    eta <- phi * gamma / rho
  }
  if(mark.influence){
    delta <- switch(case, 0.0, par.ests[4], 0.0, par.ests[5])
    xi <- par.ests[length(par.ests) - 1]
    beta <- par.ests[length(par.ests)]
    eta <- eta * (1 + delta * beta / (1 - xi))
  }
  cluster.size <- 1 / (1 - eta)
  out <- c((eta < 1), eta, cluster.size)
  names(out) <- c("stationary", "eta", "cluster.size")
  out
}
