# influence index plot  written 9 Dec 09 by S. Weisberg
# 21 Jan 10: added wrapper influenceIndexPlot(). J. Fox
# 30 March 10: bug-fixes and changed arguments, S. Weisberg
# 15 October 13:  Bug-fix on labelling x-axis

influenceIndexPlot <- function(model, ...)
	{UseMethod("infIndexPlot")}

infIndexPlot <- function(model, ...)
         {UseMethod("infIndexPlot")}

infIndexPlot.lm <- function(model,
     vars=c("Cook", "Studentized", "Bonf", "hat"), 
     main="Diagnostic Plots",
     labels, id.method = "y", 
     id.n = if(id.method[1]=="identify") Inf else 0,
     id.cex=1, id.col=palette()[1], grid=TRUE, ...) {
   what <- pmatch(tolower(vars), 
                  tolower(c("Cook", "Studentized", "Bonf", "hat")))
   if(length(what) < 1) stop("Nothing to plot")
   names <- c("Cook's distance", "Studentized residuals",
           "Bonferroni p-value", "hat-values")
# check for row.names, and use them if they are numeric.
   if(missing(labels)) labels <-  row.names(model$model)
   op <- par(mfrow=c(length(what), 1), mar=c(1, 4, 0, 2) + .0,
             mgp=c(2, 1, 0), oma=c(6, 0, 6, 0))
   oldwarn <- options()$warn
   options(warn=-1)
   xaxis <- as.numeric(row.names(model$model))
   options(warn=oldwarn)
   if (any (is.na(xaxis))) xaxis <- 1:length(xaxis)
   on.exit(par(op))
   outlier.t.test <- pmin(outlierTest(model, order=FALSE, n.max=length(xaxis),
      cutoff=length(xaxis))$bonf.p, 1)
   nplots <- length(what)
   plotnum <- 0
   for (j in what){
      plotnum <- plotnum + 1 
      y <- switch(j, cooks.distance(model), rstudent(model),
                     outlier.t.test, hatvalues(model))
      plot(xaxis, y, type="n", ylab=names[j], xlab="", xaxt="n", tck=0.1, ...)
	    if(grid){
        grid(lty=1, equilogs=FALSE)
        box()}
      if(j==3) {
            for (k in which(y < 1)) lines(c(xaxis[k], xaxis[k]), c(1, y[k]))}
          else {
            points(xaxis, y, type="h", ...)}  
      points(xaxis, y, type="p", ...)  
      if (j == 2) abline(h=0, lty=2 )
      axis(1, labels= ifelse(plotnum < nplots, FALSE, TRUE))
      showLabels(xaxis, y, labels=labels,
            id.method=id.method, id.n=id.n, id.cex=id.cex,
            id.col=id.col)
    }
    mtext(side=3, outer=TRUE ,main, cex=1.2, line=1)
    mtext(side=1, outer=TRUE, "Index", line=3)
   invisible()
}
