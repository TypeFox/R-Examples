print.eiMDsum <- function(x, digits = max(3, getOption("digits") - 4), ...) {
  cat("\nFormula: ", deparse(x$call$formula), "\n")
  cat("Total sims: ", (x$call$burnin) + (x$call$sample * x$call$thin), 
"\n")
  cat("Burnin discarded: ", x$call$burnin, "\n")
  cat("Sims saved: ", x$call$sample, "\n\n")

   "%w/o%" <- function(x,y) x[!x %in% y]
  
  if (x$short)
    cat("\nAcceptance ratios for Beta (averaged over units):\n")
  else 
    cat("\nAcceptance ratios for Beta:\n")
  print.default(format(x$acc.ratios$beta.acc, digits = digits),
                print.gap = 1, quote = FALSE)

  for (ii in names(x$acc.ratios) %w/o% c("beta.acc"))  {
    cat(paste("\nAcceptance ratios for ",
              strsplit(ii, ".acc", fixed = TRUE)[1], ":\n", sep = ""))
    print.default(format(x$acc.ratios[[ii]], digits = digits),
                  print.gap = 1, quote = FALSE)
  }
  for (ii in names(x$draws) %w/o% c("Beta", "Cell.counts")) { 
    cat(paste("\nDraws for ", ii, ":\n", sep = ""))
    print.default(format(x$draws[[ii]], digits = digits),
                  print.gap = 1, quote = FALSE)
  }
  
  cat("\nAggregate cell counts (summed over units):\n")
  print.default(format(x$draws$Cell.counts, digits = digits), print.gap = 1,
                quote = FALSE)

  invisible(x)
}
