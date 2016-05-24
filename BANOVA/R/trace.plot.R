trace.plot <-
function(x, save = FALSE){
  if (class(x) == 'BANOVA.Multinomial')
    traceplot(x$samples_l2_param, colnames(x$dMatrice$X_full[[1]]), colnames(x$dMatrice$Z), save)
  else
    traceplot(x$samples_l2_param, colnames(x$dMatrice$X), colnames(x$dMatrice$Z), save)
}
