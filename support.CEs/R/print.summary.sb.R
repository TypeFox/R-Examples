print.summary.sb <- function(x, ...) 
{
  cat("\nResults of the Street and Burgess Method\n")
  
  cat("\nOperation:\n")
  cat(x$operation, "\n")
  
  if (x$operation == "construct") {
    cat("\nTreatment combinations entered:\n")
    print(x$treatment, ...)

    cat("\nGenerators entered:\n")
    print(x$generators, ...)
  } else {
    cat("\nChoice sets entered:\n")
    print(x$choice.sets.check, ...)
  }

  cat("\nEffects entered:\n")
  cat(x$effect, "\n")

  if (x$effect != "main") {
    cat("\nDeterminant of C entered:\n")
    cat(x$determinant.c, "\n")
    if (x$effect == "mplussome") {
      cat("\nTwo-factor interactions entered:\n")
      cat(x$interaction, "\n")
    }
  }

  cat("\nMessages from processing stage:\n")
  cat(x$messages, "\n\n")

  cat(x$version, "\n")

  cat("\nChoice sets created:\n")
  print(x$choice.sets, quote = FALSE, ...)

  cat("\nB matrix:\n")
  print(x$b.chr, quote = FALSE, ...)

  cat("\nLambda matrix:\n")
  print(x$l.chr, quote = FALSE, ...)

  cat("\nC matrix:\n")
  print(x$c.chr, quote = FALSE, ...)

  cat("\nC inverse:\n")
  print(x$cinv.chr, quote = FALSE, ...)

  cat("\nCorrelation matrix:\n")
  print(x$correlation.chr, quote = FALSE, ...)

  cat("\n")

  invisible(x)
}
