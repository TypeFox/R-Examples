#-------------------------------------------------------------------------------
# Revision history:
# 2009-10-29: renamed var argument to .vcov; tidied code. John
# 2010-07-02; added method for survreg and coxph objects.
# 2010-07-02; rewrote default method to permit prarmeter names to have
#   meta-characters 
# 2011-07028  Removed meta-character checks; removed parameterPrefix because
#   it didn't work and caused problems; added parameterNames to restore the
#   utility of parameterPrefix
# 2011-10-02 Fixed bugs in the .survreg and .coxph methods so parameterNames
#   works correctly
# 2012-03-02: fixed abbreviation of envir argument. J. Fox
# 2012-04-08: modfied deltaMethod.default() to use coef and vcov
# 2012-12-10: removed the 'deltaMethodMessageFlag'
# 2013-06-20: added deltaMethod.merMod(). J. Fox
# 2013-06-20: tweaks for lme4. J. Fox
# 2013-07-01: New 'constants' argument for use when called from within a function.
# 2013-07-18: fixed a bug in passing the 'func' argument
#-------------------------------------------------------------------------------

deltaMethod <- function (object, ...) {
	UseMethod("deltaMethod")
}

deltaMethod.default <- function (object, g, vcov., func = g, constants, ...) {
	if (!is.character(g)) 
		stop("The argument 'g' must be a character string")
	if ((exists.method("coef", object, default=FALSE) ||
				(!is.atomic(object) && !is.null(object$coefficients))) 
			&& exists.method("vcov", object, default=FALSE)){
		if (missing(vcov.)) vcov. <- vcov(object)
		object <- coef(object)
	}
	para <- object         
	para.names <- names(para)
	g <- parse(text = g)
	q <- length(para)
	for (i in 1:q) {
		assign(names(para)[i], para[i])
	}
	if(!missing(constants)){
     for (i in seq_along(constants)) assign(names(constants[i]), constants[[i]])}
	est <- eval(g)
	names(est) <- NULL
	gd <- rep(0, q)
	for (i in 1:q) {
		gd[i] <- eval(D(g, names(para)[i]))
	}
	se.est <- as.vector(sqrt(t(gd) %*% vcov. %*% gd))
	data.frame(Estimate = est, SE = se.est, row.names = c(func))
}

deltaMethod.lm <- function (object, g, vcov. = vcov, 
           parameterNames = names(coef(object)), ...) {
#  if( !exists("deltaMethodMessageFlag", envir=.carEnv)){
#     message("deltaMethod arguments have changed, see help(deltaMethod)")
#     assign("deltaMethodMessageFlag", TRUE, envir=.carEnv)
#     }
	para <- coef(object)
	para.names <- parameterNames
	para.names[1] <- gsub("\\(Intercept\\)", "Intercept", para.names[1])
	names(para) <- para.names
	vcov. <- if (is.function(vcov.)) 
			vcov.(object)
		else vcov.
	deltaMethod.default(para, g, vcov.,  ...)
}

# nls has named parameters so parameterNames is ignored
deltaMethod.nls <- function(object, g, vcov.=vcov,...){
	vcov. <- if(is.function(vcov.)) vcov.(object)
	deltaMethod.default(coef(object), g, vcov., ...)   
}

deltaMethod.polr <- function(object,g,vcov.=vcov,...){
	sel <- 1:(length(coef(object)))
	vcov. <- if(is.function(vcov.)) vcov.(object)[sel, sel]
	deltaMethod.lm(object, g, vcov., ...)
}

deltaMethod.multinom <- function(object, g, vcov.=vcov, 
   parameterNames = if(is.matrix(coef(object)))
     colnames(coef(object)) else names(coef(object)), ...){
	out <- NULL
	coefs <- coef(object)
	if (!is.matrix(coefs)) { coefs <- t(as.matrix(coefs)) }
	colnames(coefs) <- parameterNames
	nc <- dim(coefs)[2]
	for (i in 1:dim(coefs)[1]){
		para <- coefs[i, ]
		ans <- deltaMethod(para, g, vcov.(object)[(i - 1) + 1:nc, (i - 1) + 1:nc], ...)
		rownames(ans)[1] <- paste(rownames(coefs)[i], rownames(ans)[1])
		out <- rbind(out,ans)
	}
	out}

# method for survreg objects. 
deltaMethod.survreg <- function(object, g, vcov. = vcov, 
           parameterNames = names(coef(object)), ...) {
 deltaMethod.lm(object, g, vcov., parameterNames , ...) }


 # method for coxph objects.
deltaMethod.coxph <- function(object, g, vcov. = vcov, 
           parameterNames = names(coef(object)), ...) {
 deltaMethod.lm(object, g, vcov.,  parameterNames, ...) }
           
# lmer

deltaMethod.merMod <- function(object, g, vcov. = vcov,
                            parameterNames = names(fixef(object)), ...) {
    deltaMethod.mer(object=object, g=g, vcov.=vcov, 
                       parameterNames=parameterNames, ...)
}
    
deltaMethod.mer <- function(object, g, vcov. = vcov,
           parameterNames = names(fixef(object)), ...) {
  para <- fixef(object)
  names(para) = parameterNames
 	vcov. <- if (is.function(vcov.)) 
			vcov.(object)
		else vcov.
  deltaMethod(para, g, vcov., ...)
  }


#lme
deltaMethod.lme <- function(object, g, vcov. = vcov,
           parameterNames = names(fixef(object)), ...) {
  para <- fixef(object)
  names(para) = parameterNames
 	vcov. <- if (is.function(vcov.)) 
			vcov.(object)
		else vcov.
  deltaMethod(para, g, vcov., ...)
  }
  
# nlsList  lsList
deltaMethod.lmList <- function(object, g, ...) {
  out <- t(sapply(object, function(x) deltaMethod(x, g, ...)))
  rownames(out) <- paste(rownames(out), g)
  out
  }
  

