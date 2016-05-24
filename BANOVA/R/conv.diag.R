conv.diag <-
function(x){
  if (class(x) == 'BANOVA.Multinomial')
    sol <- conv.geweke.heidel(x$samples_l2_param, colnames(x$dMatrice$X_full[[1]]), colnames(x$dMatrice$Z))
  else
    sol <- conv.geweke.heidel(x$samples_l2_param, colnames(x$dMatrice$X), colnames(x$dMatrice$Z))
  class(sol) <- 'conv.diag'
  sol
}
