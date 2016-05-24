### Functions for OLS and GLS Estimation
### Author: Konstantin Kashin
### August 1, 2013

ols <- function(env.base){
	coef <- coef(get("pw.output",envir=env.base)$pw.lm)
	model.mat.pw <- get("model.mat.pw",envir=env.base)
	Omega <- get("Omega",envir=env.base)
	sandwich <- t(model.mat.pw) %*% Omega %*% model.mat.pw
    vcov <- solve(t(model.mat.pw) %*% model.mat.pw) %*% sandwich %*% solve(t(model.mat.pw) %*% model.mat.pw)
    out <- list(coef=coef,vcov=vcov)
	return(out)	
}

gls <- function(env.base){
	model.mat.pw <- get("model.mat.pw",envir=env.base)
	Omega <- get("Omega",envir=env.base)
	response <- model.response(model.frame(get("pw.output",envir=env.base)$pw.lm), "numeric")    
	coef <- solve(t(model.mat.pw) %*% solve(Omega) %*% model.mat.pw) %*% t(model.mat.pw) %*% solve(Omega) %*% response
    vcov <- solve(t(model.mat.pw) %*% solve(Omega) %*% model.mat.pw)
	out <- list(coef=coef,vcov=vcov)
	return(out)	
}
