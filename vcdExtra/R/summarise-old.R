
# summarise a glm object or glmlist

summarise <- function(...) {
  .Deprecated("LRstats")
  LRstats(...)
}

# summarise <- function(object, ...) {
# 	UseMethod("summarise")
# }
# 
# stat.summarise <- function(deviance, df, onames, n) {
# 	p <- pchisq(deviance, df, lower.tail=FALSE)
# 	aic <- deviance - 2*df
# 	if (missing(n)) {
# 		result <- data.frame(aic, deviance, df, p)
# 		names(result) <- c("AIC", "LR Chisq", "Df", "Pr(>Chisq)")
# 	}
# 	else {
# 		bic <- deviance - log(n)*df
# 		result <- data.frame(aic, bic, deviance, df, p)
# 		names(result) <- c("AIC", "BIC", "LR Chisq", "Df", "Pr(>Chisq)")
# 	}
# 	
# 	rownames(result) <- onames
# 	attr(result, "heading") <- "Model Summary:"
# 	class(result) <- c("anova", "data.frame")
# 	result
# }
# 
# 
# summarise.glm <-function(object, ..., test=NULL){
# 	dotargs <- list(...)
# 	is.glm <- unlist(lapply(dotargs, function(x) inherits(x, "glm")))
# 	dotargs <- dotargs[is.glm]
# 	if (length(dotargs)) 
# 		return(summarise.glmlist(c(list(object), dotargs), test = test))
# 	
# 	oname <- as.character(sys.call())[2]
# 	result <- stat.summarise(object$deviance, object$df.residual, oname, sum(fitted(object)))
# 	result
# }
# 
# summarise.glmlist <-function(object, ..., test=NULL, sortby=NULL){
# 	nmodels <- length(object)
# 	if (nmodels == 1) 
# 		return(summarise.glm(object[[1]], test = test))
# 	if (is.null(names(object))) {
# 		oname <- as.character(sys.call())[-1]
# 		oname <- oname[1:length(object)]
# 	}
# 	else oname <- names(object)
# 	
# 	resdf <- as.numeric(lapply(object, function(x) x$df.residual))
# 	resdev <- as.numeric(lapply(object, function(x) x$deviance))
# 	n <- as.numeric(lapply(object, function(x) sum(fitted(x))))
# 	result <- stat.summarise(resdev, resdf, oname, n)
# 	if (!is.null(sortby)) {
# 		result <- result[order(result[,sortby], decreasing=TRUE),]
# 	}
# 	result
# }
# 
# 
# summarise.loglm <-function(object, ...){
# 	dotargs <- list(...)
# 	is.loglm <- unlist(lapply(dotargs, function(x) inherits(x, "loglm")))
# 	dotargs <- dotargs[is.loglm]
# 	if (length(dotargs)) 
# 		return(summarise.loglmlist(c(list(object), dotargs)))
# 	
# 	oname <- as.character(sys.call())[2]
# 	result <- stat.summarise(object$deviance, object$df, oname, sum(fitted(object)))
# 	result
# }
# 
# summarise.loglmlist <-function(object, ..., sortby=NULL){
# 	nmodels <- length(object)
# 	if (nmodels == 1) 
# 		return(summarise.loglm(object[[1]]))
# 	if (is.null(names(object))) {
# 		oname <- as.character(sys.call())[-1]
# 		oname <- oname[1:length(object)]
# 	}
# 	else oname <- names(object)
# 	
# 	resdf <- as.numeric(lapply(object, function(x) x$df))
# 	resdev <- as.numeric(lapply(object, function(x) x$deviance))
# 	n <- as.numeric(lapply(object, function(x) sum(fitted(x))))
# 	result <- stat.summarise(resdev, resdf, oname, n)
# 	if (!is.null(sortby)) {
# 		result <- result[order(result[,sortby], decreasing=TRUE),]
# 	}
# 	result
# }
# 
