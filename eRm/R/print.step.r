print.step <- function(x, ...)
{
  cat("\nResults for stepwise item elimination:\n")
  cat("Number of steps:",x$nsteps,"\n")

  if (!is.null(x$res.wald)) {
    cat("Criterion: Waldtest\n\n")
    print(round(x$res.wald, 3))
    cat("\n")
  }
  
  if (!is.null(x$res.itemfit)) {
    cat("Criterion: Itemfit\n\n")
    print(round(x$res.itemfit, 3))
    cat("\n")
  }
  
  if (!is.null(x$res.LR)) {
    cat("Criterion: Andersen's LR-test\n\n")
    print(round(x$res.LR, 3))
    cat("\n")
  }
  invisible(x)
}