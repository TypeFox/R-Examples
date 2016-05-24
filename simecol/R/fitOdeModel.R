

`fitOdeModel` <-
function(simObj, whichpar = names(parms(simObj)),
  obstime, yobs,
  sd.yobs = as.numeric(lapply(yobs, sd)),
  initialize = TRUE,
  weights = NULL,
  debuglevel = 0,
  fn = ssqOdeModel,
  method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "PORT", "newuoa", "bobyqa"),
  lower = -Inf, upper = Inf, scale.par = 1, control = list(), ...)  {

  method <- match.arg(method)

  par <- parms(simObj)[whichpar]
  if (!(min(lower) == -Inf | max(upper) == Inf)) {
    lower <- lower[whichpar]
    upper <- upper[whichpar]
  }

  upper. <-  Inf
  lower. <- -Inf
  if (!(method %in% c("L-BFGS-B", "PORT", "bobyqa"))) {
    upper. <- upper
    lower. <- lower
    upper <-   Inf
    lower <- - Inf
  }

  par <- p.unconstrain(par, lower., upper.)

  if (method == "PORT") {
    m <- nlminb(start = par, objective = fn, 
             simObj = simObj, obstime = obstime,
             yobs = yobs, sd.yobs = sd.yobs,
             pnames = names(par),  # workaround as nlminb did not pass the names  in R < 2.8.1
             initialize = initialize,
             weights = weights,
             debuglevel = debuglevel,
             scale = scale.par,
             control = control, 
             lower = lower, upper = upper)
    m$value <- m$objective # for consistency with optim
  } else if (method == "bobyqa") {
    m <- bobyqa(par = par, fn = fn, simObj = simObj, obstime = obstime,
             yobs = yobs, sd.yobs = sd.yobs,
             pnames = names(par),
             initialize = initialize,
             weights = weights,
             debuglevel = debuglevel,
             lower = lower, upper = upper,
             control = control, ...)
    names(m$par) <- names(par)  # otherwise names are not saved
    m$value      <- m$fval      # for consistency with optim
    m$message <- m$msg
  } else if (method == "newuoa") {
    m <- newuoa(par = par, fn = fn, simObj = simObj, obstime = obstime,
                yobs = yobs, sd.yobs = sd.yobs,
                pnames = names(par),
                initialize = initialize,
                lower. = lower.,
                upper. = upper.,
                weights = weights,
                debuglevel = debuglevel,
                #lower = lower, upper = upper,
                control = control, ...)
    names(m$par) <- names(par)  # otherwise names are not saved
    m$value      <- m$fval      # for consistency with optim
    m$message <- m$msg  
  } else {
    m <- optim(par, fn = fn, simObj = simObj, obstime = obstime,
             yobs = yobs, sd.yobs = sd.yobs,
             initialize = initialize,
             lower. = lower.,
             upper. = upper.,
             weights = weights,
             debuglevel = debuglevel,
             method = method,
             lower = lower, upper = upper,
             control = control, ...)
  }
  cat(m$message, "\n")
  m$par <- p.constrain(m$par, lower., upper.)
  
  obj <- new("modelFit", value=m$value, par=m$par, message=m$message, list=m)
  invisible(obj)
}

