# Modified Nov. 24, 2009 by S. Weisberg to use showLabels
#   rather than showExtremes
# 11 & 20 January 2010: changed lty=3 to lty=1 for fitted curve. J. Fox
# 14 April 2010: set id.n = 0. J. Fox
# 15 April 2010; rewrite showLabels
# 25 May 2010 added grid() to plots, S. Weisberg
# 15 August 2010, fixed so col= works correctly with plot, but not Boxplot
# 15 August 2010, deleted pch= argument, as it wasn't used
# 17 January 2011, allow spline terms; plot against
#   predict(model, type="terms")[[term.name]]
# 1 February 2011 default for AsIs changed to TRUE
# 31 March 2011 tukeyNonaddTest updated to check that yhat^2 is not 
#   a linear combination of other predictors (as in 1-way anova).
# 6 April 2011 omit printing lack-of-fit if no lack-of-fit test is possible
# 16 June 2011 allow layout=NA, in which case the layout is not set in this
#  function, so it is the responsibility of the user
# 10 Feb 2013:  adjusted colinearity check in tukeyNonaddTest
# 21 March 2013:  fixed nonconstant variance test with missing values for glms
# 11 July 2013:  wording changes
# 11 July 2013:  'groups' arg for residualPlot and residualPlots.
# 19 July 2014:  type='rstudent' fixed
# 7 October 2014: trapped error resulting from groups= when n<3

residualPlots <- function(model, ...){UseMethod("residualPlots")}

residualPlots.default <- function(model, terms= ~ . , 
     layout=NULL, ask, main="", 
     fitted=TRUE, AsIs=TRUE, plot=TRUE, tests=TRUE, groups, ...){
  mf <- if(!is.null(terms)) termsToMf(model, terms) else NULL
  groups <- if (!missing(groups)) {
      termsToMf(model, as.formula(paste("~",
           deparse(substitute(groups)))))$mf.vars[, 2, drop=FALSE]
      } else {
      if(is.null(mf$mf.groups)) NULL else
                  mf$mf.groups[, 2, drop=FALSE]
      }  
  mf <- mf$mf.vars
  vform <- update(formula(model), attr(mf, "terms"))
  if(any(is.na(match(all.vars(vform), all.vars(formula(model))))))
     stop("Only regressors in the formula can be plotted.")
  terms <- attr(mf, "term.labels") # this is a list
  vterms <- attr(terms(vform), "term.labels")
# drop interactions (order > 1)
  vterms <- setdiff(vterms, terms[attr(mf, "order") > 1])
# keep only terms that are numeric or integer or factors or poly
  good <- NULL
  for (term in vterms) if(
      (AsIs == TRUE & inherits(model$model[[term]], "AsIs")) |
      inherits(model$model[[term]], "numeric") |
      inherits(model$model[[term]], "integer") |
      (inherits(model$model[[term]], "factor") & is.null(groups)) |
      inherits(model$model[[term]], "matrix") |
      inherits(model$model[[term]], "poly")) good <- c(good, term)
  nt <- length(good) + fitted
  nr <- 0  
  if (nt == 0) stop("No plots specified")
  if (nt > 1 & plot == TRUE & (is.null(layout) || is.numeric(layout))) {
    if(is.null(layout)){
         layout <- switch(min(nt, 9), c(1, 1), c(1, 2), c(2, 2), c(2, 2), 
                             c(3, 2), c(3, 2), c(3, 3), c(3, 3), c(3, 3))
    }
    ask <- if(missing(ask) || is.null(ask)) prod(layout)<nt else ask
    op <- par(mfrow=layout, ask=ask, no.readonly=TRUE, 
            oma=c(0, 0, 1.5, 0), mar=c(5, 4, 1, 2) + .1)
    on.exit(par(op))
    }
  ans <- NULL       
  if(!is.null(good)){
    for (term in good){
     nr <- nr + 1
     qtest <- if(is.null(groups))
                 residualPlot(model, term, plot=plot, ...) else
                 residualPlot(model, term, plot=plot, groups=groups, ...)
     if(!is.null(qtest)){
        ans <- rbind(ans, qtest)
        row.names(ans)[nr] <- term}
    } }   
  # Tukey's test
  if (fitted == TRUE){      
   tuk <- if(is.null(groups))
             residualPlot(model, "fitted", plot=plot, ...) else
             residualPlot(model, "fitted", plot=plot, groups=groups, ...)
   if (!is.null(tuk)  & class(model)[1] == "lm"){
      ans <- rbind(ans, tuk)
      row.names(ans)[nr + 1] <- "Tukey test"
      ans[nr + 1, 2] <- 2*pnorm(abs(ans[nr + 1, 1]), lower.tail=FALSE)}} 
  if(plot == TRUE) mtext(side=3, outer=TRUE, main, cex=1.2)
  if(!is.null(ans)) {
     dimnames(ans)[[2]] <- c("Test stat", "Pr(>|t|)")
     return(if(tests == FALSE | !is.null(groups)) invisible(ans) else
        if(all(is.na(ans))) warning("No possible lack-of-fit tests") else 
        round(ans, 3)) } else
  invisible(NULL)
  }
  
residualPlots.lm <- function(model, ...) {
 residualPlots.default(model, ...)
 }
  
residualPlots.glm <- function(model, ...) {
 residualPlots.default(model, ...)
 }

residualPlot <- function(model, ...) UseMethod("residualPlot")

residualPlot.default <- function(model, variable = "fitted", type = "pearson",
                 groups, 
                 plot = TRUE,
                 linear = TRUE,     
                 quadratic = if(missing(groups)) TRUE else FALSE, 
                 smoother=NULL, smoother.args=list(), 
                 col.smooth=palette()[3],
                 labels, 
                 id.method = "r", 
                 id.n = if(id.method[1]=="identify") Inf else 0,
                 id.cex=1, id.col=palette()[1], 
                 col = palette()[1], col.quad = palette()[2],
                 pch=1,
                 xlab, ylab, lwd = 1, lty = 1,  
                 grid=TRUE, key=!missing(groups), ...) {
 string.capitalize <- function(string) {
     paste(toupper(substring(string, 1, 1)), substring(string, 2), sep="")}
 if(missing(labels)) 
      labels <- names(residuals(model)[!is.na(residuals(model))])
 ylab <- if(!missing(ylab)) ylab else
         paste(string.capitalize(type), "residuals")
 column <- match(variable, names(model$model))
 if(is.na(column) && variable != "fitted")
   stop(paste(variable, "is not a regressor in the mean function"))
 horiz <- if(variable == "fitted") predict(model) else model$model[[column]]
 lab <- if(variable == "fitted") {
    if(inherits(model, "glm")) 
       "Linear Predictor" else "Fitted values"} else variable
 lab <- if(!missing(xlab)) xlab else lab
 if(class(horiz)[1] == "ordered") horiz <- factor(horiz, ordered=FALSE)
 ans <-
   if(inherits(horiz, "poly")) {
       horiz <- horiz[ , 1]
       lab <- paste("Linear part of", lab)
       c(NA, NA)}
   else if (inherits(horiz, "matrix")) {
       horiz <- try(predict(model, type="terms"), silent=TRUE)
       if(class(horiz) == "try-error") 
          stop("Could not plot spline terms") 
       warning("Splines replaced by a fitted linear combination")
       horiz <- horiz[ , variable]
       c(NA, NA)
       }
   else if (inherits(horiz, "factor")) c(NA, NA)
   else residCurvTest(model, variable)
# are there groups
 if(!missing(groups)){
   if(is.data.frame(groups)){
      groups.name <- names(groups)[1]
      groups <- groups[, 1, drop=TRUE]
      }  else 
      groups.name <- deparse(substitute(groups))
   groups <- if(class(groups)[1] == "factor") groups else factor(groups, ordered=FALSE)
   if(key){ 
     mar3 <- 1.1 + length(levels(groups))
     op <- par(mar=c(5.1, 4.1, mar3, 2.1))
     on.exit(par(op))
     }   
   colors <- if(length(col) >=length(levels(groups))) col else palette()
   col <- colors[as.numeric(groups)]
   pchs <- if(length(pch) >= length(levels(groups))) pch else 1:length(levels(groups))
   pch <-  pchs[as.numeric(groups)] 
 }
 theResiduals <- switch(type, "rstudent"=rstudent(model), 
               "rstandard"=rstandard(model), residuals(model, type=type))
 if(plot==TRUE){
  if(class(horiz) == "factor") {
     idm <- if(is.list(id.method)) {
            lapply(id.method, function(x) if(x[1]=="xy") "y" else x)} else {
            if(id.method[1] == "xy") "y"}    
     Boxplot(theResiduals, horiz, xlab=lab, ylab=ylab, labels=labels, 
            id.method=idm, id.n=id.n, id.cex=id.cex,  
            id.col=id.col, ...) 
     abline(h=0, lty=2) } else 
     {    
     plot(horiz, theResiduals, xlab=lab, ylab=ylab, type="n", ...)
	   if(grid){
       grid(lty=1, equilogs=FALSE)
       box()}
     points(horiz, theResiduals, col=col, pch=pch, ...)
     if(linear){
        if(missing(groups)){abline(h=0, lty=2, lwd=2)} else {
        for (g in 1:length(levels(groups)))
             try(abline(lm(theResiduals ~ horiz, 
                       subset=groups==levels(groups)[g]), lty=2, lwd=2,
                       col=colors[g]), silent=TRUE)
        }}
     if(quadratic){
       new <- seq(min(horiz), max(horiz), length=200)
       if(missing(groups)){
          if(length(unique(horiz)) > 2){
             lm2 <- lm(theResiduals ~ poly(horiz, 2))
             lines(new, predict(lm2, list(horiz=new)), lty=1, lwd=2, col=col.quad)
             }} else {
          for (g in 1:length(levels(groups))){
             if(length(unique(horiz)) > 2){
             lm2 <- lm(theResiduals~poly(horiz, 2),
                subset=groups==levels(groups)[g])
             lines(new, predict(lm2, list(horiz=new)), lty=1, lwd=1.5, col=colors[g])
             }}}}
     if(is.function(smoother))
       if(missing(groups)){
       smoother(horiz, theResiduals, col.smooth, log.x=FALSE, log.y=FALSE,
          spread=FALSE, smoother.args=smoother.args)} else
       for (g in 1:length(levels(groups))){
          sel <- groups == levels(groups)[g]
          smoother(horiz[sel], theResiduals[sel], colors[g], log.x=FALSE, log.y=FALSE,
             spread=FALSE, smoother.args=smoother.args)}
     if(key & !missing(groups)){
       items <- paste(groups.name, levels(groups), sep= " = ")
       plotArrayLegend("top", items=items, col.items=colors, pch=pchs)
       }
     showLabels(horiz, theResiduals, labels=labels, 
            id.method=id.method, id.n=id.n, id.cex=id.cex, 
            id.col=id.col)  
        }
      }  
  invisible(ans)}
 
residCurvTest <- function(model, variable) {UseMethod("residCurvTest")}
residCurvTest.lm <- function(model, variable) {
 if(variable == "fitted") tukeyNonaddTest(model) else {
  if(is.na(match(variable, attr(model$terms, "term.labels"))))
     stop(paste(variable, "is not a term in the mean function")) else {
     xsqres <- qr.resid(model$qr, model.frame(model)[[variable]]^2)
     r <- residuals(model, type="pearson")
     m1 <- lm(r ~ xsqres, weights=weights(model))
     df.correction <- sqrt((df.residual(model)-1) / df.residual(m1))
     test <- summary(m1)$coef[2, 3] * df.correction
     c(Test=test, Pvalue=2 * pt(-abs(test), df.residual(model)-1))
     }}}

residCurvTest.glm <- function(model, variable) {
 if(variable == "fitted") c(NA, NA) else {
  if(is.na(match(variable, attr(model$terms, "term.labels"))))
     stop(paste(variable, "is not a term in the mean function")) else {
     newmod <- paste(" ~ . + I(", variable, "^2)")
     m2 <- update(model, newmod, start=NULL)
     c(Test= test<-deviance(model)-deviance(m2), Pvalue=1-pchisq(test, 1))
}}}

residCurvTest.negbin <- function(model, variable) {
  if(variable == "fitted") c(NA, NA) else {
    if(is.na(match(variable, attr(model$terms, "term.labels"))))
      stop(paste(variable, "is not a term in the mean function")) else {
        newmod <- paste(" ~ . + I(", variable, "^2)")
        m2 <- update(model, newmod, start=NULL)
        c(Test= test<-m2$twologlik - model$twologlik, Pvalue=1-pchisq(test, 1))
      }}}
     
tukeyNonaddTest <- function(model){
 tol <- model$qr$tol
 qr <- model$qr
 fitsq <- predict(model, type="response")^2
 fitsq <- qr.resid(qr, fitsq/sqrt(sum(fitsq^2)))
 if(sd(fitsq) < tol) {
    return(c(Test=NA, Pvalue=NA))
 } else {
    r <- residuals(model, type="pearson")
    m1 <- lm(r ~ fitsq, weights=weights(model))
    df.correction <- sqrt((df.residual(model) - 1)/df.residual(m1))
    tukey <- summary(m1)$coef[2, 3] * df.correction
    c(Test=tukey, Pvalue=2*pnorm(-abs(tukey)))
    }
 }
 

 
residualPlot.lm <- function(model, ...) {
  residualPlot.default(model, ...)
  }
  
residualPlot.glm <- function(model, variable = "fitted", type = "pearson", 
                 plot = TRUE, quadratic = FALSE, 
                 smoother = loessLine, smoother.args=list(k=3), ...) {
  residualPlot.default(model, variable=variable, type=type, plot=plot, 
                 quadratic=quadratic, smoother=smoother, 
                 smoother.args=smoother.args, ...)
  }

