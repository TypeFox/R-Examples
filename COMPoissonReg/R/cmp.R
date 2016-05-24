cmp <- function(formula, initial.est=NULL, nuinit=1, max=100, ...){
	initial_glm = glm(formula, family='poisson', ...)
	if (is.null(initial.est)) {
		initial.est = coef(initial_glm)
	}
	print(initial.est)

	data_frame = model.frame(formula, ...)

	response_name = names(data_frame)[1]
	x_names = names(data_frame)[-1]

	object_result = list()
	object_result$formula = formula
	object_result$response_name = response_name
	object_result$x_names = x_names
	object_result$response = data_frame[,response_name]
	object_result$predictors = as.data.frame(data_frame[,x_names])
	names(object_result$predictors) = x_names

	object_result$max = max

	internal_result = ComputeBetasAndNuHat(object_result$predictors, object_result$response, betainit=initial.est, nuinit=nuinit, max=object_result$max)
	if(internal_result$convergence == 1) {stop(sprintf("Constant CMP estimates could not be determined.  Optimization scheme did not converge: ", internal_result))}

	num_pars = length(internal_result$par)
	object_result$glm_coefficients = coef(initial_glm)
	object_result$coefficients = internal_result$par[1:(num_pars -1)]
	object_result$nu = internal_result$par[num_pars]

	attr(object_result, "class") <- c("cmp", attr(object_result, "class"))
	#class(object_result) = c("cmp", class(internal_result))
	return(object_result)
}

coef.cmp <- function(object) {
	return(object$coefficients)
}

nu.cmp <- function(object, ...) {
	return(object$nu)
}

sdev.cmp <- function(object, ...) {
	object$sdev <- CMPStdErrors(object$predictors, object$coef, object$nu, max=object$max)
	names(object$sdev) = c("(Intercept)",object$x_names, "nu")
	return(object$sdev)
}

chisq.cmp <- function(object, ...) {
	LRT(object$predictors, object$response, object$glm_coefficients, object$coef, object$nu, object$max)$teststat[1,1]
}

pval.cmp <- function(object, ...) {
	LRT(object$predictors, object$response, object$glm_coefficients, object$coef, object$nu, object$max)$pvalue
}

leverage.cmp <- function(object, ...) {
	CMPLeverage(object$predictors, object$response, object$coef, object$nu, object$max)
}

deviance.cmp <- function(object, ...) {
	CMPDeviance(object$predictors, object$response, object$coef, object$nu, leverage.cmp(object), object$max)
}

residuals.cmp <- function(object, ...) {
	return(object$response - predict(object, newdata=object$predictors))
}

predict.cmp <- function(object, ...) {
	newdata = list(...)[["newdata"]]
	return(constantCMPfitsandresids(object$coef, object$nu, newdata[,object$x_names])$fit)
}

parametric_bootstrap.cmp <- function(object, ...) {
	n = list(...)[["n"]]
	if (is.null(n)) {
		n = 1000
	}
	bootstrap_results = as.data.frame(CMPParamBoot(x=object$predictors, object$glm_coefficients, betahat=object$coef, nuhat=object$nu, n=n)$CMPresult)
	names(bootstrap_results) = c("(Intercept)",object$x_names,"nu",recursive=TRUE)
	return(bootstrap_results)
}