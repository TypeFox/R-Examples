fround <- function(x, digits)
    format(round(x, digits), nsmall=digits)

print.bugs <- function(x, digits.summary = 1, ...)
{
    if(!is.null(x$model.file))
        cat("Inference for Bugs model at \"", x$model.file, "\", ", sep="")
    if(!is.null(x$program))
        cat("fit using ", x$program, ",", sep="")
    cat("\nCurrent: ", x$n.chains, " chains, each with ", x$n.iter,
        " iterations (first ", x$n.burnin, " discarded)", sep = "")
    if(x$n.thin > 1) cat(", n.thin =", x$n.thin)
    cat("\nCumulative: n.sims =", x$n.sims, "iterations saved\n")
    print(round(x$summary, digits.summary), ...)

    if(x$n.chains > 1) {
      cat("\nFor each parameter, n.eff is a crude measure of effective sample size,")
      cat("\nand Rhat is the potential scale reduction factor (at convergence, Rhat=1).\n")
    }

    if(x$isDIC) {
      msgDICRule <- ifelse(x$DICbyR,
                           "(using the rule, pD = var(deviance)/2)", ## Gelman tweak
                           "(using the rule, pD = Dbar-Dhat)")       ## BUGS
      cat(paste("\nDIC info ", msgDICRule, "\n", sep=""))
      if(length(x$DIC) == 1) {
        cat("pD =", fround(x$pD, 1), "and DIC =", fround(x$DIC, 1))
      } else if(length(x$DIC)>1) {
        print(round(x$DIC, 1))
      }
      cat("\nDIC is an estimate of expected predictive error (lower deviance is better).\n")
    }
    invisible(x)
}
