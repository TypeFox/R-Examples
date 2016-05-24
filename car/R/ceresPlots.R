# CERES plots (J. Fox)

# last modified 9 October 2009 by J. Fox
# modified  26 Nov 2009 by S. Weisberg
#   changed layout and point marking.
#   modified 15 Mar 2010 by S. Weisberg to make the following work:
#   m1 <- lm(longley)
#   ceresPlots(longley)
# 14 April 2010: set id.n = 0. J. Fox
# new args for showLabels 15 April S. Weisberg
# modified 2 Sept 2010 by S. Weisberg, made colors, axes lables, and
# arguments more consistent with other functions; ... passes args to plot
# and boxplot.
# 16 June 2011 allow layout=NA, in which case the layout is not set in this
#  function, so it is the responsibility of the user
# 14 Sept 2012 use the ScatterplotSmoothers in car
# 18 Sept 2012 restore smooth and span args
# 20 Aug 2013 replace residuals.glm() with residuals(). John

ceresPlots<-function(model, terms= ~ ., layout=NULL, ask, main, ...){
  terms <- if(is.character(terms)) paste("~", terms) else terms
  vform <- update(formula(model), terms)
  if(any(is.na(match(all.vars(vform), all.vars(formula(model))))))
     stop("Only predictors in the formula can be plotted.")
  mf <- attr(model.frame(model), "terms")
  terms <- attr(mf, "term.labels") # this is a list
  vterms <- attr(terms(vform), "term.labels")
  good <- NULL
	if (any(attr(terms(model),"order")>1)) {
		stop("CERES plots not available for models with interactions.")}  
  for (term in vterms) if(
      inherits(model$model[[term]], "numeric") |
      inherits(model$model[[term]], "integer")) good <- c(good,term)
  nt <- length(good)
  if(length(good) < length(vterms))
    warning("Factors skipped in drawing CERES plots.")
  vterms <- good
  if (nt == 0) stop("No plots specified")
  if (missing(main)) main <- if (nt == 1) "CERES Plot" else "CERES Plots"
  if (nt == 0) stop("No plots specified")
  if (nt > 1 & (is.null(layout) || is.numeric(layout))) {
    if(is.null(layout)){
         layout <- switch(min(nt, 9), c(1, 1), c(1, 2), c(2, 2), c(2, 2), 
                             c(3, 2), c(3, 2), c(3, 3), c(3, 3), c(3, 3))
    }
    ask <- if(missing(ask) || is.null(ask)) prod(layout)<nt else ask
    op <- par(mfrow=layout, ask=ask, no.readonly=TRUE, 
            oma=c(0, 0, 1.5, 0), mar=c(5, 4, 1, 2) + .1)
    on.exit(par(op))
    }
	if(!is.null(class(model$na.action)) && 
		class(model$na.action) == 'exclude') class(model$na.action) <- 'omit'
  for(term in vterms) 
		ceresPlot(model, term, main="", ...)
	mtext(side=3, outer=TRUE, main, cex=1.2)
	invisible(0)
	}


ceresPlot<-function (model, ...) {
	UseMethod("ceresPlot")
}

ceresPlot.lm<-function(model, variable, 
  id.method = list(abs(residuals(model, type="pearson")), "x"),
  labels, 
  id.n = if(id.method[1]=="identify") Inf else 0,
  id.cex=1, id.col=palette()[1],
  line=TRUE, smoother=loessLine, smoother.args=list(), smooth, span,
	col=palette()[1], col.lines=palette()[-1],
  xlab, ylab, pch=1, lwd=2,  
  grid=TRUE, ...){
	# the lm method works with glm's too              
	if(missing(labels)) labels <- names(residuals(model))
    # smooth and span for backwards compatibility
    if (!missing(smooth)) {
        smoother <- if (isTRUE(smooth)) loessLine else FALSE
    }
    if (!missing(span)) smoother.args$span <- span
	expand.model.frame <- function (model, extras, envir = environment(formula(model)),
		na.expand = FALSE){  # modified version of R base function
		f <- formula(model)
		data <- eval(model$call$data, envir)
		ff <- foo ~ bar + baz
		if (is.call(extras)) 
			gg <- extras
		else gg <- parse(text = paste("~", paste(extras, collapse = "+")))[[1]]
		ff[[2]] <- f[[2]]
		ff[[3]][[2]] <- f[[3]]
		ff[[3]][[3]] <- gg[[2]]
		if (!na.expand) {
			naa <- model$call$na.action
			subset <- model$call$subset
 			rval <- if (is.null(data)){ 
            eval(call("model.frame", ff, data=model.frame(model),
 							   subset = subset, na.action = naa), envir)} else           
            eval(call("model.frame", ff, data = data,         
 							subset = subset, na.action = naa), envir)           
		}
		else {
			subset <- model$call$subset
			rval <- eval(call("model.frame", ff, data = data, subset = subset, 
					na.action = I), envir)
			oldmf <- model.frame(model)
			keep <- match(rownames(oldmf), rownames(rval))
			rval <- rval[keep, ]
			class(rval) <- "data.frame"
		}
		return(rval)
	}
	if(!is.null(class(model$na.action)) && 
		class(model$na.action) == 'exclude') class(model$na.action) <- 'omit'
#	if (missing(iter)){
#		iter<-if(("glm"==class(model)[1]) &&
#				("gaussian"!=as.character(family(model))[1]))
#				0
#			else 3
#	}    # use nonrobust smooth for non-gaussian glm
	var<-if (is.character(variable) & 1==length(variable)) variable
		else deparse(substitute(variable))
	mod.mat<-model.matrix(model)
	obs<-names(residuals(model))
	all.obs<-if (is.null(model$call$data)) obs else row.names(eval(model$call$data))
	xx<-rep(NA, length(all.obs))
	names(xx)<-all.obs
	terms<-predictor.names(model)
	if (is.na(match(var, terms))) stop(paste(var,"is not in the model."))
	if (!is.null(model$contrasts[[var]])) stop(paste(var,"is a factor."))
	terms<-terms[-match(var,terms)]
	if (any(attr(terms(model),"order")>1)) {
		stop("ceres plot not available for models with interactions.")
	}
	.x<-xvars<-NULL
	for (xvar in terms){
		if (is.null(model$contrasts[[xvar]])){
			xvars<-c(xvars,xvar)
			xx[obs]<-fitted.values(loess(as.formula(paste("mod.mat[,'",xvar,"']~mod.mat[,'",var,"']",sep=""))))
			.x<-cbind(.x, xx)
		}
	}
	if (is.null(xvars)) stop("There are no covariates.")
	n.x<-length(xvars)
	mf<-na.omit(expand.model.frame(model, all.vars(formula(model))))
	rownames(.x)<-all.obs
	mf$.x<-.x[obs,]
	aug.model <- update(model, . ~ . + .x, data=mf, subset=NULL)
	aug.mod.mat<-model.matrix(aug.model)
	coef<-coefficients(aug.model)
	k<-length(coef)
	posn<-k:(k-n.x+1)
	partial.res<-residuals(aug.model, "partial")[,var] +
		aug.mod.mat[,posn] %*% as.matrix(coef[posn])
	xlab <- if(!missing(xlab)) xlab else var
	ylab <- if(!missing(ylab)) ylab else
	   paste("CERES Residual(",responseName(model),")", sep="")
	plot(mod.mat[,var], partial.res, xlab=xlab, col=col, pch=pch,
		ylab=ylab, type="n", ...)
	if(grid){
    grid(lty=1, equilogs=FALSE)
    box()}
	points(mod.mat[,var], partial.res, col=col, pch=pch) 
	showLabels(mod.mat[,var], partial.res, labels=labels, 
            id.method=id.method, id.n=id.n, id.cex=id.cex,
            id.col=id.col)
	if (line) abline(lm(partial.res~mod.mat[,var]), lty=2, lwd=lwd, 
            col=col.lines[1])
	if (is.function(smoother)) {
    smoother(mod.mat[, var], partial.res, col=col.lines[2], log.x=FALSE,
       log.y=FALSE, spread=FALSE, smoother.args=smoother.args)
#		lines(lowess(mod.mat[,var], partial.res, iter=iter, f=span), lwd=lwd,
#          col=col.lines[2])
	}
}                    

ceresPlot.glm<-function(model, ...){
	ceresPlot.lm(model, ...)
}
