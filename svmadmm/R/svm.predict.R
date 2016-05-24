#' The predict function for svm.admm
#' 
#' @description
#'   \code{svm.predict} The function applies a model produced by the \code{svm.admm} function
#'    to every row of a data matrix and returns the model predictions.
#' @param x.te A n*p data matrix, input test data.
#' @param model The training result created by \code{svm.admm}.
#' @return n-length vector, predicted labels.
#' 
#' @export 
svm.predict <- function(x.te, model){
	# if (!requireNamespace("kernlab")) {
	# 	install.packages("kernlab")
	# 	# library("kernlab")
 #  		# stop('The package kernlab was not installed')
	# }
	n = nrow(x.te)
	if (model$type == 0)
	{
		rbf1 = model $ kern
		k.te = kernlab::kernelMatrix(rbf1, x.te, model$x.tr)
		fit = as.vector(sign(sum(model$alpha * model$y.tr) + k.te %*% (model$y.tr * model$alpha)))
	}
	else
	{
		fit = as.vector(sign(cbind(rep(1, n), x.te) %*% model$beta))
	}
	return(fit)
}
