print.boot.lmf <-
function(x,
                           digits = max(3, getOption("digits") - 3),
                           ...)
{
  #Title
  cat("\nESTIMATING FLUCTUATING SELECTION IN AGE-STRUCTURED POPULATIONS\n",
      sep = "")
  #Display numbers of bootstraps generated
  cat(sprintf(ngettext(x$nboot, "%s Bootstrap replicate generated\n",
                       "%s Bootstrap replicates generated\n"), paste(x$nboot)), sep = "")
  #Display the call
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  #Print running time and optim running time
  cat("Running time:\n")
  print.default(format(c(total = x$running.time, optim = x$optim.time),
                       digits = digits), print.gap = 2, quote = FALSE)
  #End of print
  cat("-End-\n")
  #Return 'x' invisibly
  invisible(x)
}
