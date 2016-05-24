# generic functions defined in rms
#
#   contrast <- function(fit, ...) UseMethod("contrast")
#   gendata <- function(fit, ...) UseMethod("gendata")

# this overrides the version in nlme which doesn't work for us
formula.gls <- function(x, ..., env=NULL)
{
   if (is.null(env))
      # this seems silly, but it's the way gls defines formula
      eval(x$call$model)
   else
      # this gives you some control if you know what you're doing
      eval(x$call$model, env)
}

# methods neeeded, but missing, from geepack
coef.geese <- function(object, ...) object$beta
vcov.geese <- function(object, ...) object$vbeta

# This function mimics rms:::gendata, which is unusable
# on non-rms objects
generateData <- function(fit, factors, ..., env=NULL)
{
   tt <- tryCatch(terms(fit), error=function(e) terms(formula(fit, env=env)))
   order <- attr(tt, 'order')
   tlabs <- attr(tt, 'term.labels')
   nam <- tlabs[order == 1]
   fnam <- names(factors)
   nf <- length(factors)

   if (nf == 0)
      stop('illegal factors argument')

   wh <- charmatch(fnam, nam, 0)
   if (any(wh == 0))
      stop(paste("factor(s) not in design:", paste(names(factors)[wh == 0], collapse=" ")))

   if (nf < length(nam))
      stop('not enough factors')

   settings <- list()
   for (i in 1:nf)
      settings[[fnam[i]]] <- factors[[i]]

   expand.grid(settings)
}





