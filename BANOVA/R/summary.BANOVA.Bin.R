summary.BANOVA.Bin <-
function(object, ...){
  sol <- list(anova.table = object$anova.table,
              coef.table = object$coef.tables$full_table,
              pvalue.table = object$pvalue.table, conv = object$conv, call = object$call)
  class(sol) <- "summary.BANOVA.Bin"
  sol
}
