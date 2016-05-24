# October 1, 2014 mcPlots, by S. Weisberg and J. Fox
# 'mc' stands for Marginal and Conditional:  for the specified regressor X in model 
# The 'marginal' plot is of Y vs X with Y and X both centered
# The 'conditional plot is the added-variable plot e(Y|rest) vs e(X|rest)
# If 'overlaid=TRUE', the default, both plots are overlayed
# If 'overlaid=FALSE', then the plots are side-by-side
# The 'overlaid' plot is similar to the initial and final frame of an ARES plot
# Cook and Weisberg (1989), "Regression diagnostics with dynamic graphics", Technometrics, 31, 277.
# This plot would benefit from animation.

mcPlots <- function(model, terms=~., layout=NULL, ask, overlaid=TRUE, ...){
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
#	if (intercept) good <- c("(Intercept)", good)
  if(attr(attr(model.frame(model), "terms"), "intercept") == 0)
    stop("Error---the 'lm' object must have an intercept")
	nt <- length(good)
	if (nt == 0) stop("No plots specified")
#	if (missing(main)) main <- if (nt == 1) paste("Marginal/Conditional Plot:", good) else 
#                                                "Marginal/Conditional Plots"
  if (nt == 0) stop("No plots specified")
  if(overlaid){
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
  } else{
    if (nt >= 1 & (is.null(layout) || is.numeric(layout))) {
      if(is.null(layout)){
        layout <- switch(min(nt, 4), c(1, 2), c(2, 2), c(3, 2), c(4, 2))
      }
      ask <- if(missing(ask) || is.null(ask)) layout[1] < nt else ask
      op <- par(mfrow=layout, ask=ask, no.readonly=TRUE, 
                oma=c(0, 0, 1.5, 0), mar=c(5, 4, 1, 2) + .1)
      on.exit(par(op))
    }
  }
	for (term in good) mcPlot(model, term, new=FALSE, overlaid=overlaid, ...)
#	mtext(side=3,outer=TRUE,main, cex=1.2)
}



mcPlot <-  function(model, ...) UseMethod("mcPlot")

mcPlot.lm <- function(model, variable,
                      id.method = list(abs(residuals(model, type="pearson")), "x"),
                      labels, 
                      id.n = if(id.method[1]=="identify") Inf else 0,
                      id.cex=1, id.col=palette()[1],
                      col.marginal="blue", col.conditional="red", col.arrows="gray",
                      pch = c(16,1), lwd = 2, grid=TRUE,   ###removed arg main
                      ellipse=FALSE, ellipse.args=list(levels=0.5), 
                      overlaid=TRUE, new=TRUE, ...){
  variable <- if (is.character(variable) & 1 == length(variable))
                 variable  else deparse(substitute(variable))
  if(new && !overlaid) {
    op <- par(mfrow=c(1,2))
    on.exit(par(op))
  }
  if(missing(labels)) 
    labels <- names(residuals(model)[!is.na(residuals(model))])
  else deparse(substitute(variable))
  if(attr(attr(model.frame(model), "terms"), "intercept") == 0)
    stop("Error---the 'lm' object must have an intercept")
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
  res0 <- lm(cbind(mod.mat[, var], response) ~ 1, weights=wt)$residual
  res  <- lsfit(mod.mat[, -var], cbind(mod.mat[, var], response),
               wt = wt, intercept = FALSE)$residuals
  xlab <- paste(var.names[var], "| others") 
  ylab <- paste(responseName, " | others")  
  xlm <- c( min(res0[, 1], res[, 1]), max(res0[, 1], res[, 1]))
  ylm <- c( min(res0[, 2], res[, 2]), max(res0[, 2], res[, 2]))
  if(overlaid){ 
     plot(res[, 1], res[, 2], xlab = xlab, ylab = ylab, type="n", 
          main=paste("Marginal/Conditional plot of", var.names[var]),
          xlim=xlm, ylim=ylm,  ...)
     if(grid){
       grid(lty=1, equilogs=FALSE)
       box()}     
     points(res0[, 1], res0[, 2], pch=pch[1], col=col.marginal)
     points(res[, 1], res[, 2], col=col.conditional, pch=pch[2], ...)
     arrows(res0[, 1], res0[, 2], res[, 1], res[, 2], length=0.125, col=col.arrows)
     abline(lsfit(res0[, 1], res0[, 2], wt = wt), col = col.marginal, lwd = lwd)
     abline(lsfit(res[, 1], res[, 2], wt = wt), col = col.conditional, lwd = lwd)
     if (ellipse) {
       ellipse.args1 <- c(list(res0, add=TRUE, plot.points=FALSE, col=col.marginal), ellipse.args)
       do.call(dataEllipse, ellipse.args1)
       ellipse.args1 <- c(list(res, add=TRUE, plot.points=FALSE, col=col.conditional), ellipse.args)
       do.call(dataEllipse, ellipse.args1)
     }
     showLabels(res0[, 1],res0[, 2], labels=labels, 
              id.method=id.method, id.n=id.n, id.cex=id.cex, 
              id.col=id.col)
    colnames(res) <- c(var.names[var], responseName)
    rownames(res) <- rownames(mod.mat)
    invisible(res)} 
  else { # side.by.side plots
    plot(res0[, 1], res0[, 2], type="n", 
          xlab = paste("Centered", var.names[var], sep=" "), 
          ylab = paste("Centered", responseName, sep=" "), 
          main=paste("Marginal plot of", var.names[var]),
          xlim=xlm, ylim=ylm,  ...)
    if(grid){
       grid(lty=1, equilogs=FALSE)
       box()}     
    points(res0[, 1], res0[, 2], pch=pch[1], col=col.marginal)
    abline(lsfit(res0[, 1], res0[, 2], wt = wt), col = col.marginal, lwd = lwd)
    if (ellipse) {
      ellipse.args1 <- c(list(res0, add=TRUE, plot.points=FALSE, col=col.marginal), ellipse.args)
      do.call(dataEllipse, ellipse.args1)
    }
    showLabels(res0[, 1],res0[, 2], labels=labels, 
               id.method=id.method, id.n=id.n, id.cex=id.cex, 
               id.col=id.col)
    colnames(res) <- c(var.names[var], responseName)
    rownames(res) <- rownames(mod.mat)    
    plot(res[, 1], res[, 2], xlab = xlab, ylab = ylab, type="n", 
         main=paste("Added-Variable plot of", var.names[var]),
         xlim=xlm, ylim=ylm,  ...)
    if(grid){
      grid(lty=1, equilogs=FALSE)
      box()}
    points(res[, 1], res[, 2], col=col.conditional, pch=pch[2], ...)
    abline(lsfit(res[, 1], res[, 2], wt = wt), col = col.conditional, lwd = lwd)
    if (ellipse) {
       ellipse.args1 <- c(list(res, add=TRUE, plot.points=FALSE, col=col.conditional), ellipse.args)
       do.call(dataEllipse, ellipse.args1)
     }
    showLabels(res[, 1],res[, 2], labels=labels, 
                id.method=id.method, id.n=id.n, id.cex=id.cex, 
                id.col=id.col)
    invisible(res)}      
}

