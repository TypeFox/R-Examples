make_cvdglars <- function(object,control){
	beta <- object$b
	names(beta) <- c("Int.",colnames(object$X))
	var_cv <- names(beta[abs(beta) > 0])[-1]
	out <- list(call=NULL,family=NULL,formula_cv=NULL,var_cv=var_cv,beta=beta,dev_m=object$dev_m,dev_v=object$dev_v,g0=control$g0,
				g_hat=object$g0,g_max=object$g[1],X=object$X,y=object$y,conv=object$conv,control=control)
	class(out) <- "cvdglars"
	out
}
