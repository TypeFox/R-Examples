# Bootstrap sampling
lexpit.bootstrap <- function(object, nboot = 100){
	
bootstrap.sampling <- function(nboot = 50, data, case=NULL, strata=NULL){
	
	sampling <- function(data, case=NULL, strata=NULL){
		index <- 1:nrow(data)
	unlist( tapply(index, list(case, strata), function(x) sample(x, length(x), replace = TRUE)) )
	}

sapply(1:nboot, function(x) sampling(data=data, case=case, strata=strata))
}

bootstrap.fit <- function(object, index){

	fit <- tryCatch(lexpit(object@formula.linear, object@formula.expit,
			  data = object@data[index,], 
			  weights = object@weights[index],
			  strata = object@strata[index]),
			  error = function(x) NA)

	if(class(fit)!="lexpit")
		rep(NA, length(coef(object)))
	else	
		coef(fit)
}


	samples <- bootstrap.sampling(nboot=nboot, data=object@data, case = object@y, strata = object@strata)
	result <- apply(samples, 2, function(x) bootstrap.fit(object = object, index = x))
	
	keep <- !is.na(result[1,])
	result <- result[,keep]
	coefs <- coef(object)
	n <- object@p+object@q
	
	Q <- apply(result, 2, function(x) outer(x-coefs, x-coefs))
	Q <- rowMeans(Q)
	Q <- matrix(Q, n, n)
	
	object@vcov.linear <- matrix(Q[1:object@p,1:object@p])
	object@vcov.expit <- matrix(Q[(object@p+1):n,(object@p+1):n])

object
}

