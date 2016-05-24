SUE.fitted.values <-
function(fit){
mf=model.frame(formula=fit$p$formula, data=fit$p$data)
z=(as.matrix(mf[,-1])%*%as.matrix(fit$coeff[-1])+fit$coeff[1])
z
}
