print.stratasize <-
function(x,...)
{
  cat("\nstratamean object: Stratified sample size determination\n", sep="")
  cat("\ntype of sample: ", x$call$type, "\n", sep="")
  cat("\ntotal sample size determinated: ", x$n, "\n", sep="")
  invisible(x)
}
