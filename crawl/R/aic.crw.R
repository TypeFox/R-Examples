#' Calculates AIC for all objects of class crwFit listed as arguments
#' 
#' AIC, delta AIC, and Akaike weights for all models listed as arguments.
#' 
#' 
#' The function can either be executed with a series of 'crwFit' objects (see
#' \code{\link{crwMLE}}) without the '.crwFit' suffix or the function can be
#' called without any arguments and it will search out all 'crwFit' objects in
#' the current workspace and produce the model selection table for all 'crwFit'
#' objects in the workspace. Caution should be used when executing the function
#' in this way. ALL 'crwFit' objects will be included whether ot not the same
#' locations are used!  For all of the models listed as arguments (or in the
#' workspace), AIC, delta AIC, and Akaike weights will be calculated.
#' 
#' @param \dots a series of crwFit objects
#' @return A table, sorted from lowest AIC value to highest.
#' @author Devin S. Johnson
#' @export
"aic.crw" <- function(...)
{
  lnms <- NULL
  models <- list(...)
  if(length(models) == 0) {
      lnms <- list()
      lx <- ls(envir=parent.frame(2))
      for (i in 1:length(lx)) {
          classval <- class(eval(parse(text=lx[i]), envir=parent.frame(2)))
          if("crwFit" %in% classval) lnms <- append(lnms,list(lx[i]))
      }
      models <- eval(parse(text=paste("list(",
                             paste(paste(lnms, "=", lnms, sep=""),
                                   collapse=","), ")")), envir=parent.frame())
      vnms <- do.call("c", lnms)
  } else {
      models <- list(...)
      vnms <- all.vars(match.call())
  }
  num.mod <- length(vnms)
  AIC.vec <- numeric(num.mod)
  ks <- numeric(num.mod)
  for (i in 1:num.mod) {
      AIC.vec[i] <- round(models[[i]]$aic, 2)
      ks[i] <- length(models[[i]]$fixPar) - sum(!is.na(models[[i]]$fixPar))
  }
  deltaAIC <- round(AIC.vec - min(AIC.vec), 2)
  wAIC <- round(exp(-0.5 * deltaAIC) / sum(exp(-0.5 * deltaAIC)), 2)
  ord <- order(deltaAIC)
  out <- data.frame(Name=vnms, k=ks, AIC=AIC.vec,
                    dAIC=deltaAIC, weight=wAIC)
  return(out[ord, ])
}




