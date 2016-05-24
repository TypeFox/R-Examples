enet.selection.A.b <-
function(y, X, lambda, alpha=1, mu){
  n = length(y)
  enet.fit = fit.enet.fixed.lambda(y=y, X=X, lambda=lambda*n, alpha=1, mu=mu)
  Ab.obj = enet.selection.A.b.from.enet(enet.fit)
  
  list(which.col=Ab.obj$which.col, A=Ab.obj$A, b=Ab.obj$b)
}
