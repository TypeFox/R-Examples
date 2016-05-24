summary.BANOVA.ordMultinomial <-
function(object, ...){
  res <- list(anova.table = object$anova.table,
              coef.table = rbind(object$coef.tables$full_table, object$coef.tables$cutp_table),
              pvalue.table = object$pvalue.table, conv = object$conv, call = object$call)
  class(res) <- "summary.BANOVA.ordMultinomial"
  res
}
