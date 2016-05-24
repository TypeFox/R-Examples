## Methods for standard generic functions (defined in package 'stats')
## for objects for which they are missing from original packages.

##------------------------------------------------------------------------------
# default 'family' method
##------------------------------------------------------------------------------
`family.default` <-
function (object, ...)  {
	cl <- get_call(object)
	if(is.null(cl)) 
		return(NULL)
	fam <- cl$family
	if(is.null(fam)) 
		fam <- formals(match.fun(cl[[1L]]))$family
	if(is.null(fam))
		return(NA)
	switch(mode(fam), call = eval(fam), name =, character = match.fun(fam)())
}

##------------------------------------------------------------------------------
## package 'glmmML'
##------------------------------------------------------------------------------

# this replaces the original method, merely to get rid of the annoying behaviour
# in summary.glmML. it does not do anything except for printing the model
# output.
`summary.glmmML` <-
function(object, ...) object

##------------------------------------------------------------------------------
## package 'nlme'
##------------------------------------------------------------------------------

`family.gls` <-
`family.lme` <-
function (object, ...) {
	if (inherits(object$family, "family")) object$family else gaussian()
}


`model.frame.lme` <-
function (formula, random = FALSE, ...) {
	x <- formula
	frm <- formula(x)
	if(random) {
		for(reStruct in x$modelStruct$reStruct)
			frm[[3L]] <- call("+", frm[[3L]], attr(reStruct, "formula")[[2L]])
	}
	mfArgs <- list(formula = frm, data = x$data[rownames(x$fitted), ], drop.unused.levels = TRUE)
	do.call("model.frame", mfArgs)
	#droplevels(do.call("model.frame", mfArgs))
}

`model.matrix.lme` <-
function (object, random = FALSE, ...) {
	mf <- model.frame(object, random = random)
	model.matrix(formula(terms(mf)), mf, contrasts.arg = object$contrasts)
}

##------------------------------------------------------------------------------
## package 'betareg'
##------------------------------------------------------------------------------

family.betareg <-
function (object, ...) {
	ret <- binomial(object$link$mean)
	attr(ret, "link-precision") <- object$link$precision
	ret
}

##------------------------------------------------------------------------------
## package 'coxme'
##------------------------------------------------------------------------------

`formula.coxme` <-
function(x, ...)  {
	ret <- x$formulaList$fixed
	f <- ret[[3L]]
	for(f1 in x$formulaList$random) f <- call("+", f, f1)
	ret[[3L]] <- f
	ret
}

`formula.lmekin` <-
function(x, ...) eval.parent(x$call$formula)

##------------------------------------------------------------------------------
## package 'pscl'
##------------------------------------------------------------------------------

`family.zeroinfl` <-
function(object, ...) binomial(link = object$link)

##------------------------------------------------------------------------------
## package 'aod'
##------------------------------------------------------------------------------

`formula.glimML` <-
function(x, ...) x@formula

`family.glimML` <-
function(object, ...) switch(object@method,
	"BB" = binomial(object@link),
	#"NB" = MASS::negative.binomial(theta = 1/object@param['phi.(Intercept)'],
	"NB" = get("negative.binomial", asNamespace("MASS"))(
		theta = 1 / object@param['phi.(Intercept)'], link = object@link))

`terms.glimML` <-
function (x, ...) terms.formula(x@formula, ...)

`model.frame.glimML` <-
function (formula, ...)
	model.frame(formula@formula, data = formula@data, na.action = formula@na.action)

##------------------------------------------------------------------------------
## package 'aod3'
##------------------------------------------------------------------------------

model.matrix.aodml <-
function (object, ...) object$X.b

`model.frame.aodml` <-
function (formula, ...)
model.frame(formula$formula, data = formula$dat)

##------------------------------------------------------------------------------
## package 'unmarked'
##------------------------------------------------------------------------------

#setMethod("logLik", "unmarkedFit", logLik.unmarkedFit)

`formula.unmarkedFit` <- function (x, ...) x@formula


##------------------------------------------------------------------------------
## package 'geepack'
##------------------------------------------------------------------------------

`coef.geese` <- 
function (object, ...) object$beta

## What if 'data' changed in the meantime?
# model.matrix.gee <-
# function (object, ...) {
	# cl <- get_call(fgee)
	# cl[[1L]] <- as.name("model.matrix")
	# cl$object <- cl$formula
	# cl$id <- cl$corstr <- cl$formula <- NULL
	# eval.parent(cl)
# }

##------------------------------------------------------------------------------
## package 'yags'
##------------------------------------------------------------------------------

`coef.yagsResult` <-
function (object, ...)
structure(object@coefficients, names = object@varnames)


`getCall.yagsResult` <-
	function(x, ...) x@Call
	
`formula.yagsResult` <-
function (x, ...) 
eval.parent(x@Call$formula)

##------------------------------------------------------------------------------
## package 'MCMCglmm'
##------------------------------------------------------------------------------

`formula.MCMCglmm` <-
function (x, ...) x$Fixed$formula


`family.MCMCglmm` <-
function (object, ...) object$family

##------------------------------------------------------------------------------
## package 'caper'
##------------------------------------------------------------------------------

`formula.caic` <-
function(x, ...) formula(x$mod)

##------------------------------------------------------------------------------
## package 'asreml'
##------------------------------------------------------------------------------

#XXX: this is for fixed effects only (should sparse be included too?)
`formula.asreml` <- 
function (x, ...)  as.formula(x$fixed.formula)


`family.asreml` <- 
function(object, ...) {
	fam <- object$family
	fam$linkfun <- fam$link
	fam$link <- fam$family[2L]
	fam$family <- fam$family[1L]
	fam$linkinv <- fam$inverse
	fam$inverse <- NULL
	class(fam) <- "family"
	fam
}

##------------------------------------------------------------------------------
## package 'maxlike'
##------------------------------------------------------------------------------

formula.maxlikeFit <-
function (x, ...) 
as.formula(get_call(x)$formula, env = parent.frame())

##------------------------------------------------------------------------------
## package 'geeM'
##------------------------------------------------------------------------------

model.matrix.geem <-
function (object, ...) 
object$X

## EOF
