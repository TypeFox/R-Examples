`DIC` <-
function (object, ...) {
	if (!missing(...)) {
		lls <- sapply(list(object, ...), function(x) {
			c(extractDIC(x), attr(logLik(x), "df"))
		})
		val <- data.frame(df = lls[2L, ], DIC = lls[1L, ])
		Call <- match.call()
		row.names(val) <- make.unique(as.character(Call[-1L]))
		val
	} else extractDIC(object)
}

if(!exists("extractDIC", mode = "function")) {
	extractDIC <- function (fit, ...) UseMethod("extractDIC")
}

## from package 'arm'
`extractDIC.merMod` <- function (fit, ...) {
	dev <- deviance(fit, REML = isREML(fit))
    devML <- deviance(fit, REML = FALSE)
    as.vector(2 * devML - dev)
}

`extractDIC.MCMCglmm` <- function (fit, ...) fit$DIC

`extractDIC.lme` <- function (fit, ...) {
	ll <- as.vector(logLik(fit, REML = fit$method == "REML"))
    llML <- as.vector(logLik(fit, REML = FALSE))
    2 * ll - 4 * llML
}








