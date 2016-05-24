
# Utility functions (J. Fox)

# 16 March 2010 changed 'vars' argument to 'terms'
# 28 June 2010 added df.terms.surveg and model.matrix.survreg
# 15 November 2010 added squeezeBlanks
# 21 January 2011 added functions to support mixed models
# 2012-04-08 added exists.method
# 2012-06-23: added call to globalVariables(). John
# 2012-12-10: added .carEnv to avoid warnings in R > 2.16.0
# 2013-06020: added .merMod methods to df.residual() and has.intercept(). John
# 2014-05-16: added .multinom method for has.intercept(). John
# 2014-08-19: added package.installed() function, unexported. John
# 2014-11-02: termsToMf fixed, Sandy
# 2015-01-13: fixed model.matrix.lme() to work with model with formula as object. John
# 2015-01-27: .carEnv now lives in the global environment. John
# 2015-09-04: added model.matrix.coxme() and alias.coxme(). John
# 2015-09-11: added some support for VGAM::vglm objects. John

#if (getRversion() >= "2.15.1") globalVariables(c(".boot.sample", ".boot.indices"))

.carEnv <- new.env(parent=globalenv())

# function to find "nice" numbers

nice <- function(x, direction=c("round", "down", "up"), lead.digits=1){
	direction <- match.arg(direction)
	if (length(x) > 1) return(sapply(x, nice, direction=direction, lead.digits=lead.digits))
	if (x == 0) return(0)
	power.10 <- floor(log(abs(x),10))
	if (lead.digits > 1) power.10 <- power.10 - lead.digits + 1
	lead.digit <- switch(direction,
		round=round(abs(x)/10^power.10),
		down=floor(abs(x)/10^power.10),
		up=ceiling(abs(x)/10^power.10))
	sign(x)*lead.digit*10^power.10
}

has.intercept <- function (model, ...) {
	UseMethod("has.intercept")
}

has.intercept.default <- function(model, ...) any(names(coefficients(model))=="(Intercept)")

has.intercept.multinom <- function(model, ...) {
  nms <- names(coef(model))
  any(grepl("\\(Intercept\\)", nms))
}
  
term.names <- function (model, ...) {
	UseMethod("term.names")
}

term.names.default <- function (model, ...) {
	term.names <- labels(terms(model))
	if (has.intercept(model)) c("(Intercept)", term.names)
	else term.names
}

predictor.names <- function(model, ...) {
	UseMethod("predictor.names")
}

predictor.names.default <- function(model, ...){
	predictors <- attr(terms(model), "variables")
	as.character(predictors[3:length(predictors)])
}

responseName <- function (model, ...) {
	UseMethod("responseName")
}

responseName.default <- function (model, ...) deparse(attr(terms(model), "variables")[[2]])

response <- function(model, ...) {
	UseMethod("response")
}

response.default <- function (model, ...) model.response(model.frame(model))

is.aliased <- function(model){
	!is.null(alias(model)$Complete)
}

df.terms <- function(model, term, ...){
	UseMethod("df.terms")
}

df.terms.default <- function(model, term, ...){
	if (is.aliased(model)) stop("Model has aliased term(s); df ambiguous.")
	if (!missing(term) && 1 == length(term)){
		assign <- attr(model.matrix(model), "assign")
		which.term <- which(term == labels(terms(model)))
		if (0 == length(which.term)) stop(paste(term, "is not in the model."))
		sum(assign == which.term)
	}
	else {
		terms <- if (missing(term)) labels(terms(model)) else term
		result <- numeric(0)
		for (term in terms) result <- c(result, Recall(model, term))
		names(result) <- terms
		result
	}
}

df.terms.multinom <- function (model, term, ...){
	nlev <- length(model$lev)
	if (!missing(term) && 1 == length(term)) {
		assign <- attr(model.matrix(model), "assign")
		which.term <- which(term == labels(terms(model)))
		if (0 == length(which.term))
			stop(paste(term, "is not in the model."))
		sum(assign == which.term) * (nlev - 1)
	}
	else {
		terms <- if (missing(term))
				labels(terms(model))
			else term
		result <- numeric(0)
		for (term in terms) result <- c(result, Recall(model,
					term))
		names(result) <- terms
		result
	}
}

df.terms.polr <- function (model, term, ...){
	if (!missing(term) && 1 == length(term)) {
		assign <- attr(model.matrix(model), "assign")
		which.term <- which(term == labels(terms(model)))
		if (0 == length(which.term))
			stop(paste(term, "is not in the model."))
		sum(assign == which.term)
	}
	else {
		terms <- if (missing(term))
				labels(terms(model))
			else term
		result <- numeric(0)
		for (term in terms) result <- c(result, Recall(model,
					term))
		names(result) <- terms
		result
	}
}

df.terms.survreg <- function(model, term, ...){
	if (is.aliased(model)) stop("Model has aliased term(s); df ambiguous.")
	if (!missing(term) && 1 == length(term)){
		assign <- attr(model.matrix(model, data=model.frame(model)), "assign")
		which.term <- which(term == labels(terms(model)))
		if (0 == length(which.term)) stop(paste(term, "is not in the model."))
		sum(assign == which.term)
	}
	else {
		terms <- if (missing(term)) labels(terms(model)) else term
		result <- numeric(0)
		for (term in terms) result <- c(result, Recall(model, term))
		names(result) <- terms
		result
	}
}

model.matrix.survreg <- function(object, ...) model.matrix.default(object, model.frame(object))

mfrow <- function(n, max.plots=0){
	# number of rows and columns for array of n plots
	if (max.plots != 0 && n > max.plots)
		stop(paste("number of plots =",n," exceeds maximum =", max.plots))
	rows <- round(sqrt(n))
	cols <- ceiling(n/rows)
	c(rows, cols)
}

inv <- function(x) solve(x)

coefnames2bs <- function(g, para.names, parameterPrefix="b"){
	metas <- c("(", ")", "[", "]", "{", "}", ".", "*", "+", "^", "$", ":", "|")
	metas2 <- paste("\\", metas, sep="")
	metas3 <- paste("\\\\", metas, sep="")
	for (i in seq(along=metas))
		para.names <- gsub(metas2[i], metas3[i], para.names) # fix up metacharacters
	para.order <- order(nchar(para.names), decreasing=TRUE) 
	para.names <- para.names[para.order] # avoid partial-name substitution
	std.names <- if ("(Intercept)" %in% para.names)
				paste(parameterPrefix, 0:(length(para.names) - 1), sep = "")
			else paste(parameterPrefix, 1:length(para.names), sep = "")
	std.names.ordered <- std.names[para.order]
	for (i in seq(along=para.names)){
		g <- gsub(para.names[i], std.names.ordered[i], g) 
	}
	list(g=g, std.names=std.names)
}


showLabelsScatter <- function(x, y, labels, id.var = NULL,
	id.method = c("mahal", "identify", "none"), log="", id.cex=.75, id.n=3, id.col=palette()[1], 
	range.x=range(.x), show=TRUE) {
	id.method <- match.arg(id.method)
	if (id.method == "none" || id.n == 0 || !show) return(invisible(NULL))
	if(id.n > 0L) {
		if (missing(labels))
			labels <- if (!is.null(id.var)) names(id.var)
				else as.character(seq(along=x))
		getPoints <- function(z) {
			names(z) <- labels
			iid <- seq(length=id.n)
			zs <- z[order(-z)[iid]]
			match(names(zs), labels)
		}
		logged <- function(axis=c("x", "y")){
			axis <- match.arg(axis)
			0 != length(grep(axis, log))
		}
		valid <- complete.cases(x, y)
		x <- x[valid]
		y <- y[valid]
		labels <- labels[valid]
		if (length(id.var) == length(valid)) 
			id.var <- id.var[valid]
		.x <- if (logged("x")) log(x) else x
		.y <- if (logged("y")) log(y) else y
		ind <- if (!is.null(id.var)) {
				if (length(id.var) == length(x)) order(-abs(id.var))[1L:id.n] 
				else if(is.character(id.var)) match(id.var, labels) else id.var
			}
			else switch(id.method,
					x = getPoints(abs(.x - mean(.x))),
					y = getPoints(abs(.y - mean(.y))),
					xy = union(getPoints(abs(.x - mean(.x))),
						getPoints(abs(.y - mean(.y)))),
					mahal= getPoints(rowSums(qr.Q(qr(cbind(1, .x, .y))) ^ 2)))
		ind <- na.omit(ind)
		if (length(ind) == 0) return(invisible(NULL))
		labpos <- c(4, 2)[1 + as.numeric(.x > mean(range.x))]
		text(x[ind], y[ind], labels[ind], cex = id.cex, xpd = TRUE,
			pos = labpos[ind], offset = 0.25, col=id.col)
		return(labels[ind])
	} 
}


#  outerLegend, written by S. Weisberg Feb 2010
#  outerLegend function
#  puts a legend in the margin, either at the upper left (margin = 3)
#  the default or upper right side otherwise
#  all the args from legend are used except for x, y, and xpd which are
#  set in the function.
#  offset is a fraction of the plot width or height to locate the legend
outerLegend <- function(..., margin=3, offset=0, adjust=FALSE){
   lims <- par("usr")
   if (margin == 3) {
      x0 <- lims[1] + offset*(lims[2]-lims[1])
      y0 <- lims[4] }
   else {
      x0 <- lims[2] + offset*(lims[2]-lims[1])
      y0 <- lims[4]
   }
   leg <- legend(x0, y0, ... , xpd=TRUE, plot=FALSE)
   if (margin == 3) {
      y0 <- y0 + leg$rect$h
      if(adjust == TRUE) x0 <- x0 - leg$text$x[1]
   }
   legend(x0, y0, ... , xpd=TRUE)
   }
   
# added by J. Fox 18 Nov 2010

squeezeBlanks <- function(text){
	gsub(" *", "",  text)
}

# added by J. Fox 21 Jan 2011 to support mixed models

df.residual.mer <- function(object, ...) NULL

df.residual.merMod <- function(object, ...) NULL

df.residual.lme <- function(object, ...) Inf

has.intercept.mer <- function(model){
	any(names(fixef(model))=="(Intercept)")
}

has.intercept.merMod <- function(model){
    any(names(fixef(model))=="(Intercept)")
}
	
model.matrix.lme <- function(object, ...){
	model.matrix(formula(object), eval(object$call$data))
}

# added by J. Fox 2012-04-08 to use in deltaMethod.default()

exists.method <- function(generic, object, default=TRUE, strict=FALSE){
	classes <- class(object)
	if (default) classes <- c(classes, "default")
	if (strict) classes <- classes[1]
	any(paste(generic, ".", classes, sep="") %in%
					as.character(methods(generic)))
}

     
# Used by marginalModelPlots, residualPlots added 2012-09-24
plotArrayLegend <- function(
      location=c("top", "none", "separate"),
      items, col.items, lty.items, lwd.items, title="legend",
      pch=1:length(items)) {
   if(location== "none") return()
   n <- length(items)
   if(location == "top" ) { # add legend
      usr <- par("usr")
      coords <-list(x=usr[1], y=usr[3])
      leg <- legend( coords, items,
                col=col.items, pch=pch,
                bty="n", cex=1, xpd=NA, plot=FALSE)
      coords <- list(x = usr[1], y=usr[4] + leg$rect$h)
      legend( coords, items,
         col=col.items, pch=pch, bty="n", cex=1, xpd=NA)
  }
  if(location == "separate") {
    plot(0:1, 0:1, xaxt="n", yaxt="n", xlab="", ylab="", type="n")
    bg <- par()$bg
    legend("center", items,
          lty=lty.items, lwd=lwd.items, fill=col.items, border=col.items,,
          col=col.items, box.col=par()$bg,
          title=title)
  }
}

termsToMf <- function(model, terms){ 
  gform <- function(formula) {
    if (is.null(formula)){
        return(list(vars=formula, groups=NULL))
    }
    has.response <- length(formula) == 3
    rhs <- if(has.response) formula[[3]] else formula[[2]]
# either a single variable or '.' on the RHS
    if (length(rhs) == 1){
        return(list(vars=formula, groups=NULL))
    }
    if (length(rhs) != 3) stop("incorrectly formatted 'terms' argument")
# check for "|", indicating grouping
    if(deparse(rhs[[1]] == "|")) {
      if(length(rhs[[3]]) > 1) stop("Only one conditioning variable is permitted")
      groups <- as.formula(paste("~", deparse(rhs[[3]])))
      rhs <- rhs[[2]]
    } else groups <- NULL 
    vars <- as.formula(paste("~ ", deparse(rhs)))
    list(vars=vars, groups=groups)
  }
  terms <- gform(as.formula(terms))
  mf.vars <- try(update(model, terms$vars, method="model.frame"),
     silent=TRUE)
# This second test is used for models like m1 <- lm(longley) which
# fail the first test because update doesn't work
  if(class(mf.vars) == "try-error")
       mf.vars <- try(update(model, terms$vars,
               method="model.frame", data=model.frame(model)), silent=TRUE)
  if(class(mf.vars) == "try-error") stop("argument 'terms' not interpretable.")
  if(!is.null(terms$groups)){
     mf.groups <- try(update(model, terms$groups, method="model.frame"), silent=TRUE)
     if(class(mf.groups) == "try-error")
       mf.groups <- try(update(model, terms$groups,
               method="model.frame", data=model.frame(model)), silent=TRUE)
     if(class(mf.groups) == "try-error") stop("argument 'terms' not interpretable.")
  } else {mf.groups <- NULL}
  list(mf.vars=mf.vars, mf.groups=mf.groups)
  }

# the following function isn't exported, tests for existance of a package:

package.installed <- function(package){
  package <- as.character(substitute(package))
  result <- try(find.package(package), silent=TRUE)
  !class(result) ==  "try-error"
}

# support for coxme objects

model.matrix.coxme <- function(object, ...){
    if (!requireNamespace("survival")) stop("survival package is missing")
    class(object) <- "coxph"
    model.matrix(object)
}

alias.coxme <- function(model){
    if(any(which <- is.na(coef(model)))) return(list(Complete=which))
    else list()
}

# to make linearHypothesis() work again and to make Anova() work with VGAM:"vglm" objects 

# df.residual.vglm <- function(object, ...) object@df.residual

# vcov.vglm <- function(object, ...) vcovvlm(object, ...)

# coef.vglm <- function(object, ...) coefvlm(object, ...)

has.intercept.vlm <- function(model, ...) any(grepl("^\\(Intercept\\)", names(coef(model))))

# formula.vglm <- function(x, ...) formulavlm(x = x, ...)

# model.matrix.vglm <- function(object, ...) model.matrixvlm(object, ...)

