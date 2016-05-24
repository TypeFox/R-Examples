smac = function(x, y, loss = c("logistic", "psvm", "boost"), 
	weight=NULL, nlambda = 100, lambda.min = ifelse(nobs < np, 0.05, 1e-03), 
	lambda = NULL, standardize = TRUE, epsilon = 1e-05)
{
	if (!is.matrix(x) & !is.data.frame(x)) stop("The covariates should be either a matrix or a data.frame.")	
	if (is.na(sum(x)) | is.nan(sum(x))) stop("There should be no NA/NaN in the covariates.")	

	np=ncol(x)
	nobs=nrow(x)

	if (length(y)!=nobs) {stop("The dimension of covariates should match the length of the label.")}
	if (is.na(epsilon) | is.nan(epsilon)) stop("Epsilon should not be NA/NaN.")
	if (!is.numeric(epsilon)) {stop("Epsilon should be numeric.")}
	if (epsilon<=0) {stop("Epsilon should be strictly positive.")}
	if (length(as.vector(epsilon))>1) {stop("Epsilon should be a scalar.")}

	if (is.na(lambda.min) | is.nan(lambda.min)) stop("lambda.min should not be NA/NaN.")
	if (!is.numeric(lambda.min)) {stop("lambda.min should be numeric.")}
	if (lambda.min>=1) {stop("lambda.min should be strictly less than 1.")}
	if (length(as.vector(lambda.min))>1) {stop("lambda.min should be a scalar.")}
	
	if (is.na(nlambda) | is.nan(nlambda)) stop("nlambda should not be NA/NaN.")
	if (!is.numeric(nlambda)) {stop("nlambda should be numeric.")}
	if (nlambda!=round(nlambda)) {stop("nlambda should be an integer.")}
	if (nlambda<1) {stop("nlambda should be at least 1.")}
	if (length(as.vector(nlambda))>1) {stop("nlambda should be a scalar.")}
	
	if (is.null(weight)) {weight = rep(1,nobs)}
	if (!is.numeric(weight)) {stop("The weight vector should be numeric.")}
	if (is.na(sum(weight)) | is.nan(sum(weight))) stop("There should be no NA/NaN in the weight vector.")
	if (min(weight)<0) {stop("The weight vector should be non-negative.")}
	if (length(weight)!=nobs) {stop("The length of the weight should agree with the number of observations.")}

	loss = match.arg(loss)
	this.call = match.call()

	fit = switch(loss, 
	logistic=angle_logi(x,y,weight,nlambda,lambda.min,lambda,standardize,epsilon),
	psvm=angle_psvm(x,y,weight,nlambda,lambda.min,lambda,standardize,epsilon),
	boost=angle_boost(x,y,weight,nlambda,lambda.min,lambda,standardize,epsilon)
			)

	fit$call = this.call
	class(fit) = "smac" 
	return(fit)
}
