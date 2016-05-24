`nobs.glmmML` <- function(object, ...) length(object$coefficients) +
	object$cluster.null.df

# Extends: nlme --- No longer needed
`nobs.gls` <- function(object, nall = TRUE, ...) {
	p <- object$dims$p
	N <- object$dims$N
	if (nall) return (N)
	REML <- object$method == "REML"
	N - REML * p
}

`nobs.lme` <- function(object, nall = TRUE, ...) {
	N <- object$dims$N
	if (nall) return (N)
	p <- object$dims$ncol[object$dims$Q + 1L]
	REML <- object$method == "REML"
	N - REML * p
}
# # p - the number of coefficients in the linear model.
# #N - the number of observations in the data,
# #Q - the number of grouping levels
# #ncol - the number of columns in the model matrix for each level of grouping from innermost to outermost
# #  (last two values are equal to the number of fixed effects and one).

# Extends: survival
`nobs.survreg` <-
function (object, ...)
length(object$linear)


# Extends: nnet/spdep
`nobs.sarlm` <-
`nobs.spautolm` <-
`nobs.multinom` <-
function(object, ...) NROW(fitted(object))

`nobs.rq` <-
function (object, ...) length(object$y)

`nobs.coxme` <-
function (object, ...) object$n[2L]

`nobs.lmekin` <-
function (object, ...) object$n[1L]

`nobs.hurdle` <-
`nobs.zeroinfl` <- `nobs.lmekin`

`nobs.glimML` <- 
function (object, ...) attr(logLik(object), "nobs")

`nobs.unmarkedFit` <- 
function(object, ...) 
	#get("sampleSize", asNamespace("unmarked"))(object)
	unmarked::sampleSize(object)

`nobs.yagsResult` <-
function (object, ...) length(object@residuals)

`nobs.asreml` <- 
`nobs.aodml` <-
`nobs.splm` <- 
function (object, ...) length(resid(object))

`nobs.MCMCglmm` <-
function (object, ...) object$Residual$nrl


`nobs.gamm` <-
	function (object, ...) getFrom("stats", "nobs.glm")(object$gam, ...)

`nobs.mark` <- 
function (object, ...) object$results[['n']]

`nobs.coxph` <- 
`nobs.pgls` <-
`nobs.logistf` <-
`nobs.phylolm` <- 
function (object, ...)
object$n

`nobs.caic` <-
function (object, ...)
nobs(object$mod)

`nobs.aodql` <-
function (object, ...) 
nobs(object$fm)

`nobs.cplm` <-
function (object, ...) 
sum(!is.na(resid(object)))

`nobs.cpglmm` <-
function (object, ...) 
object@dims[['n']]

`nobs.maxlikeFit` <-
function (object, ...)
nrow(object[['points.retained']])

`nobs.geem` <-
function (object, ...) 
sum(object$weights != 0)
