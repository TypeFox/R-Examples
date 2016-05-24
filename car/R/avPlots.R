# October 23, 2009  avPlots by S. Weisberg.  avPlot by John Fox
# 13 January 2010: changed default id.n=3. J. Fox
# 13 March 2010: added intercept argument. J. Fox
# 14 April 2010: set id.n = 0. J. Fox
# 22 April 2010: modified id.n S. Weisberg
# 10 May 2010:  added gridlines
# 25 May 2010:  changed default color scheme
# 5 June 2011: made several modifications, slightly adapting code contributed by M. Friendly (J. Fox):
#   added ellipse, ellipse.args arguments 
#   added main argument to avPlot.lm
#   return x and y residuals invisibly
# 16 June 2011 allow layout=NA, in which case the layout is not set in this
#  function, so it is the responsibility of the user
# 22 Sept 2013 added argument marginal.scale to set xlim and ylim according to xlim and 
#  ylim of marginal plot (S. Weisberg)



avPlots <- function(model, terms=~., intercept=FALSE, layout=NULL, ask, 
		main, ...){
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
	if (missing(main)) main <- if (nt == 1) paste("Added-Variable Plot:", good) else "Added-Variable Plots"
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
	res <- as.list(NULL)
	for (term in good) res[[term]] <- avPlot(model, term, main="", ...)
	mtext(side=3,outer=TRUE,main, cex=1.2)
	invisible(res)
}

avp <- function(...) avPlots(...)

avPlot <-  function(model, ...) UseMethod("avPlot")

avPlot.lm <- function(model, variable,
		id.method = list(abs(residuals(model, type="pearson")), "x"),
		labels, 
		id.n = if(id.method[1]=="identify") Inf else 0,
		id.cex=1, id.col=palette()[1],
		col = palette()[1], col.lines = palette()[2],
		xlab, ylab, pch = 1, lwd = 2, main=paste("Added-Variable Plot:", variable), grid=TRUE,
		ellipse=FALSE, ellipse.args=NULL, 
    marginal.scale=FALSE, ...){
	variable <- if (is.character(variable) & 1 == length(variable))
				variable
			else deparse(substitute(variable))
	if(missing(labels)) 
		labels <- names(residuals(model)[!is.na(residuals(model))])
	else deparse(substitute(variable))
	mod.mat <- model.matrix(model)
	var.names <- colnames(mod.mat)
	var <- which(variable == var.names)
	if (0 == length(var))
		stop(paste(variable, "is not a column of the model matrix."))
	response <- response(model)
	responseName <- responseName(model)
	if (is.null(weights(model)))
		wt <- rep(1, length(response))
	else wt <- weights(model)
	res <- lsfit(mod.mat[, -var], cbind(mod.mat[, var], response),
			wt = wt, intercept = FALSE)$residuals
	xlab <- if(missing(xlab)) paste(var.names[var], "| others") else xlab
	ylab <- if(missing(ylab)) paste(responseName, " | others")  else ylab
  adjrange <- function(z, zres, wt){
    wtmean <- sum(z*wt)/sum(wt)
    c(min(z-wtmean, zres), max(z-wtmean, zres))
  }
  xlm <- if(marginal.scale) adjrange(mod.mat[, var], res[, 1], wt) else NA
	ylm <- if(marginal.scale) adjrange(response, res[, 2], wt) else NA
  if(marginal.scale) {
    plot(res[, 1], res[, 2], xlab = xlab, ylab = ylab, type="n", main=main, xlim=xlm, ylim=ylm, ...) 
  } else{
	  plot(res[, 1], res[, 2], xlab = xlab, ylab = ylab, type="n", main=main, ...)
  }
	if(grid){
		grid(lty=1, equilogs=FALSE)
		box()}
	points(res[, 1], res[, 2], col=col, pch=pch, ...)
	abline(lsfit(res[, 1], res[, 2], wt = wt), col = col.lines, lwd = lwd)
	if (ellipse) {
		ellipse.args <- c(list(res, add=TRUE, plot.points=FALSE), ellipse.args)
		do.call(dataEllipse, ellipse.args)
	}
	showLabels(res[, 1],res[, 2], labels=labels, 
			id.method=id.method, id.n=id.n, id.cex=id.cex, 
			id.col=id.col)
	colnames(res) <- c(var.names[var], responseName)
	rownames(res) <- rownames(mod.mat)
	invisible(res)  
}

avPlot.glm <- function(model, variable, 
		id.method = list(abs(residuals(model, type="pearson")), "x"),
		labels,
		id.n = if(id.method[1]=="identify") Inf else 0,
		id.cex=1, id.col=palette()[1], 
		col = palette()[1], col.lines = palette()[2],
		xlab, ylab, pch = 1, lwd = 2,  type=c("Wang", "Weisberg"), 
		main=paste("Added-Variable Plot:", variable), grid=TRUE,
		ellipse=FALSE, ellipse.args=NULL, ...) {
	type<-match.arg(type)
	if(missing(labels)) labels <- names(residuals(model)[!is.na(residuals(model))])
	variable<-if (is.character(variable) & 1==length(variable)) variable
			else deparse(substitute(variable))
	mod.mat<-model.matrix(model)
	var.names<-colnames(mod.mat)
	var<-which(variable==var.names)
	if (0==length(var)) stop(paste(variable,"is not a column of the model matrix."))
	response<-response(model)
	responseName<-responseName(model)
	wt<-model$prior.weights
	mod<-glm(response~mod.mat[,-var]-1, weights=wt, family=family(model))
	res.y<-residuals(mod, type="pearson")
	wt<-if (type=="Wang") wt*model$weights else wt
	res.x<-lsfit(mod.mat[,-var], mod.mat[,var], wt=wt,    
			intercept=FALSE)$residuals
	xlab <- if(missing(xlab)) paste(var.names[var], "| others") else xlab
	ylab <- if(missing(ylab)) paste(responseName, " | others") else ylab
	plot(res.x, res.y, xlab=xlab, type="n", ylab=ylab, main=main, ...)
	if(grid){
		grid(lty=1, equilogs=FALSE)
		box()}
	points(res.x, res.y, col=col, pch=pch, ...)
	abline(lsfit(res.x, res.y, wt=wt), col=col.lines, lwd=lwd)
	if (ellipse) {
		ellipse.args <- c(list(res.x, res.y, add=TRUE, plot.points=FALSE), ellipse.args)
		do.call(dataEllipse, ellipse.args)
	}
	showLabels(res.x,res.y, labels=labels, 
			id.method=id.method, id.n=id.n, id.cex=id.cex, 
			id.col=id.col)
	res <- cbind(res.x, res.y)
	colnames(res) <- c(var.names[var], responseName)
	rownames(res) <- rownames(mod.mat)
	invisible(res)
}

