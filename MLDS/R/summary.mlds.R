`summary.mlds` <-
function(object, 
			digits = max(3, getOption("digits") - 4), ...) {
#object, obj of class mlds
	z <- object
	cat
	ans <- list()
	ans$pscale <- z$pscale
	names(ans$pscale) <- z$stimulus
	ans$sigma <- z$sigma
	ans$logLik <- if (z$method %in% c("optim", "formula" )) 
		z$logLik else
		logLik(z$obj)[1]
	ans$method <- z$method
	if(z$method == "formula") {
		ans$formula <- z$formula
		ans$pars = z$par}
	ans$link <- z$link
	class(ans) <- "summary.mlds"
	ans
	}

`summary.mlbs` <-
function(object, 
			digits = max(3, getOption("digits") - 4), ...) {
#object, obj of class mlds
	z <- object
	cat
	ans <- list()
	ans$pscale <- z$pscale
	names(ans$pscale) <- z$stimulus
	ans$sigma <- z$sigma
	ans$logLik <- if (z$method == "glm") logLik(z$obj)[1] else logLik(z)
	ans$method <- z$method
	if(z$method == "formula") {
		ans$formula <- z$formula
		ans$pars = z$par}
	ans$link <- z$link
	class(ans) <- "summary.mlbs"
	ans
	}