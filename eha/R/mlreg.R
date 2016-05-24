mlreg <-
function (formula = formula(data),
          data = parent.frame(),
          na.action = getOption("na.action"),
          init = NULL,
          method = c("ML", "MPPL"),
          control = list(eps = 1e-8,
          maxiter = 10, n.points = 12, trace = FALSE),
          singular.ok = TRUE,
          model = FALSE,
          center = TRUE,
          x = FALSE,
          y = TRUE,
          boot = FALSE,
          geometric = FALSE,
          rs = NULL,
          frailty = NULL,
          max.survs = NULL)
{
    return("'mlreg' is deprecated; use 'coxreg' instead (see 'methods')")
    if (method[1] == "ML") method <- "ml"
    else if (method[1] == "MPPL") method <- "mppl"
    else stop(paste("Unknown method", as.character(method[1])))

    efrac <- 0
    coxreg(formula,
           data,
           t.offset = NULL,
           weights = NULL,
           na.action,
           init,
           method,
           control,
           singular.ok,
           model,
           center,
           x,
           y,
           boot,
           efrac,
           geometric,
           rs,
           frailty,
           max.survs)
}
