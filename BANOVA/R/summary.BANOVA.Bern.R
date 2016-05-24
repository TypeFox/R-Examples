summary.BANOVA.Bern <-
function(object, ...){
  res <- list(anova.table = object$anova.table,
              coef.table = object$coef.tables$full_table,
              pvalue.table = object$pvalue.table, conv = object$conv, call = object$call)
  class(res) <- "summary.BANOVA.Bern"
  res
}
