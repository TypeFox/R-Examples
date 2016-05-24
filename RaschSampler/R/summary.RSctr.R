"summary.RSctr" <-
function(object,...)
{
  cat("\nCurrent sampler control specifications in ",deparse(substitute(object)),":\n", sep="")
  cat("\tburn_in =",object$burn_in,"\n")
  cat("\tn_eff =",object$n_eff,"\n")
  cat("\tstep =",object$step,"\n")
  cat("\tseed =",object$seed,"\n")
  cat("\ttfixed =",object$tfixed,"\n\n")
  invisible(object)
}

