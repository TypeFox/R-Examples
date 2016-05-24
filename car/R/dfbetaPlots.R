# added 13 March 2010 by J. Fox
# modified 2 Sept 2010 by J. Fox, made colors, axes lables, and
# arguments more consistent with other functions; ... passes args to plot

dfbetasPlots <- function(model, ...){
	UseMethod("dfbetasPlots")
}

dfbetasPlots.lm <- function(model, terms= ~ ., intercept=FALSE, layout=NULL, ask, 
		main, xlab, ylab, labels=rownames(dfbeta), 
		id.method="y",  
		id.n=if(id.method[1]=="identify") Inf else 0, id.cex=1, 
		id.col=palette()[1], col=palette()[1], grid=TRUE, ...){
	terms <- if(is.character(terms)) paste("~",terms) else terms
	vform <- update(formula(model),terms)
	if(any(is.na(match(all.vars(vform), all.vars(formula(model))))))
		stop("Only predictors in the formula can be plotted.")
	terms.model <- attr(attr(model.frame(model), "terms"), "term.labels")
	terms.vform <- attr(terms(vform), "term.labels")
	terms.used <- match(terms.vform, terms.model)
	mm <- model.matrix(model) 
	model.names <- attributes(mm)$dimnames[[2]]
	model.assign <- attributes(mm)$assign
	good <- model.names[!is.na(match(model.assign, terms.used))]
	if (intercept) good <- c("(Intercept)", good)
	nt <- length(good)
	if (nt == 0) stop("No plots specified")
	if (missing(main)) main <- if (nt == 1) "dfbetas Plot" else "dfbetas Plots"
	if (missing(xlab)) xlab <- "Index"
	autolabel <- missing(ylab)
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
	dfbetas <- dfbetas(model)
	for (term in good) {
		dfbs <- dfbetas[, term]
		if (autolabel) ylab=term
		plot(dfbs, xlab=xlab, ylab=ylab, type="n", ...)
		if(grid){
			grid(lty=1, equilogs=FALSE)
			box()}
		points(dfbs,  col=col, ...)
		abline(h=c(-1, 0, 1), lty=2)
		showLabels(seq(along=dfbs), dfbs, id.method=id.method, 
				id.n=id.n, labels=labels, id.col=id.col,
				id.cex=id.cex, ...)
	}
	mtext(side=3,outer=TRUE,main, cex=1.2)
	invisible(NULL)
}

dfbetaPlots <- function(model, ...){
	UseMethod("dfbetaPlots")
}

dfbetaPlots.lm <- function(model, terms=~., intercept=FALSE, layout=NULL, ask, 
		main, xlab, ylab,
		labels=rownames(dfbeta), id.method="y",
		id.n=if(id.method[1]=="identify") Inf else 0, id.cex=1, 
		id.col=palette()[1], col = palette()[1], grid=TRUE, ...){
	terms <- if(is.character(terms)) paste("~",terms) else terms
	vform <- update(formula(model),terms)
	if(any(is.na(match(all.vars(vform), all.vars(formula(model))))))
		stop("Only predictors in the formula can be plotted.")
	terms.model <- attr(attr(model.frame(model), "terms"), "term.labels")
	terms.vform <- attr(terms(vform), "term.labels")
	terms.used <- match(terms.vform, terms.model)
	mm <- model.matrix(model) 
	model.names <- attributes(mm)$dimnames[[2]]
	model.assign <- attributes(mm)$assign
	good <- model.names[!is.na(match(model.assign, terms.used))]
	if (intercept) good <- c("(Intercept)", good)
	nt <- length(good)
	if (nt == 0) stop("No plots specified")
	if (missing(main)) main <- if (nt == 1) "dfbeta Plot" else "dfbeta Plots"
	if (missing(xlab)) xlab <- "Index"
	autolabel <- missing(ylab)
  if (nt > 1 &  (is.null(layout) || is.numeric(layout))) {
    if(is.null(layout)){
         layout <- switch(min(nt, 9), c(1, 1), c(1, 2), c(2, 2), c(2, 2), 
                             c(3, 2), c(3, 2), c(3, 3), c(3, 3), c(3, 3))
    }
    ask <- if(missing(ask) || is.null(ask)) prod(layout)<nt else ask
    op <- par(mfrow=layout, ask=ask, no.readonly=TRUE, 
            oma=c(0, 0, 1.5, 0), mar=c(5, 4, 1, 2) + .1)
    on.exit(par(op))
    }
	dfbeta <- dfbeta(model)
	seb <- sqrt(diag(vcov(model)))
	for (term in good) {
		dfb <- dfbeta[, term]
		se <- seb[term]
		if (autolabel) ylab <- term
		plot(dfb, xlab=xlab, ylab=ylab, type="n", ...)
		if(grid) grid(lty=1, equilogs=FALSE) 
		points(dfb, col=col, ...)
		abline(h=c(-se, 0, se), lty=2)
		showLabels(seq(along=dfb), dfb, id.method=id.method, id.n=id.n, 
				labels=labels, id.cex=id.cex, id.col=id.col, ...)
	}
	mtext(side=3,outer=TRUE,main, cex=1.2)
	invisible(NULL)
}
