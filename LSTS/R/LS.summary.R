LS.summary = function(object){
  aux1 = object$coef
  aux2 = sqrt(diag(object$var.coef))
  aux3 = aux1/aux2
  aux4 = 2*(1-pnorm(abs(aux3)))
  Table = cbind(aux1, aux2, aux3, aux4)
  colnames(Table) = c("Estimate", "Std. Error", "z-value", "Pr(>|z|)")
  list(summary = round(Table,4), aic = object$aic, npar = length(aux1))
}