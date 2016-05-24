`coeffs` <-
function (model) UseMethod("coeffs")

`coeffs.gls` <-
function (model) summary(model)$coefficients

`coeffs.lme` <-
function(model) model$coefficients$fixed

`coeffs.merMod` <-
function (model) lme4::fixef(model)

`coeffs.coxme` <-
`coeffs.lmekin` <-
function(model) {
	# for class coxme:
	ret <- model$coefficients
	# for class lmekin and older coxme
	if(is.list(ret) && !is.null(ret$fixed)) return(ret$fixed)
	ret
}

`coeffs.unmarkedFit` <- 
function(model) {
	ret <- lapply(model@estimates@estimates, coef, altNames = FALSE)
	pfx <- rep(vapply(model@estimates@estimates, slot, "", "short.name"),
		vapply(ret, length, 1L))
	ret <- unlist(unname(ret))
	Ints <- which(names(ret) == "Int")
	names(ret) <- paste0(pfx, "(", names(ret), ")")
	attr(ret, "Intercept") <- Ints
	ret
}

`coeffs.splm` <- 
function (model) {
	c(model$coefficients, model$arcoef,
	if(is.matrix(model$errcomp)) model$errcomp[, 1L] else model$errcomp)
}

`coeffs.MCMCglmm` <-
function (model)
#summary(model)$solutions[, 1L]
colMeans(model$Sol[, seq.int(model$Fixed$nfl), drop = FALSE])

`coeffs.gamm` <-
function (model) coef(model$gam)

`coeffs.mark` <- function(model) {
	cf <- model$results$beta[, 1L]
	names(cf) <- gsub("^([a-zA-Z]+):(.*)$", "\\1(\\2)",
		rownames(model$results$beta), perl = TRUE)
	cf
}

`coeffs.multinom` <- 
function (model) {
	cf <- coef(model)
	if (!is.vector(cf)) {
		cf <- t(as.matrix(cf))
    	cfnames <- expand.grid(dimnames(cf), stringsAsFactors = FALSE)
		cfnames <- sprintf("%s(%s)", cfnames[,2L], cfnames[,1L])
		structure(as.vector(cf), names = cfnames)
	} else cf
}

`coeffs.asreml` <- 
function (model) {
	coef(model)$fixed  ## should include also '$sparse' ?
}

`coeffs.cpglmm` <-
function (model) 
model@fixef


`coeffs.survreg` <-
function (model) {
	  rval <- coef(model)
	  if (nrow(vcov(model)) > length(rval)) { # scale was estimated
			lgsc <- log(model$scale)
			names(lgsc) <- if(is.null(names(lgsc)))
				  "Log(scale)" else			  
				  paste0("Log(scale):", names(lgsc))
			rval <- c(rval, lgsc)
	  }
	  rval
}

`coeffs.default` <-
#function(model) (if(isS4(model)) coef else coef)(model)
function(model) coef(model)
