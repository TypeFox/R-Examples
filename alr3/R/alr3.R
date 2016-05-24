#  R version of the alr3 package  -- a few functions not moved to car

###########################################################################
# Chapter 5  pureErrorAnova
###########################################################################
# completely rewritten 2/24/2006
# renamed pureErrorAnova
pureErrorAnova <- function(mod){UseMethod("pureErrorAnova")}
pureErrorAnova.lm <- function(mod) {
 if (is.null(mod$model)) mod <- update(mod, model=TRUE)
 p <- dim(mod$model)[2] -1
 mod$model$Lack.of.Fit <- 
   factor(randomLinComb(model.matrix(mod), 101319853))
 aov1 <- anova(mod)
 #set.seed(save.seed) # restore random number seed
 if (length(levels(mod$model$Lack.of.Fit)) == length(mod$model$Lack.of.Fit)) 
  aov1 else {
  aov2 <- anova(lm(mod$model[ , 1]~mod$model$Lack.of.Fit, weights=weights(mod)))
  rrow <- dim(aov1)[1]
  aov2[1, 1] <- aov1[rrow, 1]-aov2[2, 1]
  aov2[1, 2] <- aov1[rrow, 2]-aov2[2, 2]
  aov2[1, 3] <- aov2[1, 2]/aov2[1, 1]
  aov1[1:(rrow-1), 4] <- aov1[1:(rrow-1), 3]/aov2[2, 3]
  aov2[1, 4] <- aov2[1, 3]/aov2[2, 3]
  row.names(aov2) <- c(" Lack of fit", " Pure Error")
  aov <- rbind(aov1, aov2)
  aov[ , 5] <- pf(aov[ , 4], aov[ , 1], aov2[2, 1], lower.tail=FALSE)
  aov
  }}

 
randomLinComb <- function(X, seed=NULL) {UseMethod("randomLinComb")}

randomLinComb.default <- function(X, seed=NULL) {
 if(!is.null(seed)) set.seed(seed)
 std <- function(x){ 
    s <- sd(x)
    if( s > 0) (x-mean(x))/s else x}
 as.vector(apply(X, 2, std)%*% as.vector(2*rnorm(dim(X)[2])-1) )
 }
 
randomLinComb.lm <- function(X, ...) {
  randomLinComb(model.matrix(X), ...)}
 
randomLinComb.lm <- function(X, seed=NULL) {
 if(is.null(X$model)) X <- update(X, model=TRUE)
 randomLinComb(X$model[ , -1], seed=seed)}


###########################################################################
# Chapter 6  Delta Method
###########################################################################
# in file deltaMethod.R in 'car'

##################################################################
# pod models; Cook and Weisberg (2004), American Statistician
##################################################################

# This code was written by Lexin Li.  It was modified and simplified 
# by S. Weisberg to work with standard R stuff in nls.

pod <- function(x, ...) {UseMethod("pod")}

pod.lm <- function(x, group, mean.function, control, ...) {
 call <- x$call
 call[[1]] <- as.name("pod")
 if(!missing(group)) call$group <- as.name(as.character(substitute(group)))
 if(!missing(mean.function)) call$mean.function <- 
     as.character(substitute(mean.function))
 if(!missing(control)) call$control <- 
     as.name(as.character(substitute(control)))
 eval(call)}

pod.formula <-
function (formula, data=sys.parent(), group, subset, weights=NULL, na.action, 
    mean.function=c("pod", "common", "parallel", "general"), 
    singular.ok = FALSE, contrasts = NULL, offset, control=nls.control(), ...) 
{ 
    op <- options(warn=-1)  # suppress warnings
    mod <- match.arg(mean.function)
    call <- match.call(expand.dots=FALSE)
    call$... <- NULL
    subset <- if (missing(subset)) NULL else call$subset
    g <- factor(eval(substitute(group), data))
    gname <- substitute(group)  
    gname <- as.character(if(length(gname) == 1) gname else gname[3])
    if (mod != "pod") { 
      if (mod == "common") {
       call$group <- call$mean.function <- NULL
       call[[1]] <- as.name("lm")
       ans <- eval(call, parent.frame())} else
     if (mod == "parallel") {
       assign(gname, g)
       call$formula <- update(as.formula(call$formula), 
                       as.formula(paste("~.+", gname, sep="")))
                        #as.formula("~.+g"))
       call$group <- call$mean.function <- NULL
       call[[1]] <- as.name("lm")
       ans <- eval(call)} else
     if (mod == "general") {
       assign(gname, g)
       call$formula <- update(as.formula(call$formula), 
                        as.formula(paste("~(.)*", gname, sep="")))
       call$group <- call$mean.function <- NULL
       call[[1]] <- as.name("lm")
       ans <- eval(call)}
     ans$pod.mean.function <- mod
     ans$group <- structure(if (!is.null(subset)) g[eval(subset)] else g,  
                            name=gname)
     class(ans) <- c("pod.lm", "lm")
     ans} else
    { # POD model
    l<-nlevels(g) 
    cl <- call
    cl$group <- cl$mean.function <- NULL
    cl[[1]] <- as.name("lm") 
    fit1 <- eval(cl)
    if(fit1$rank != length(fit1$coef))
     stop("Predictors are linearly related.  Try again after deleting predictors.")
    coef1 <- coef(fit1)
    p1 <- predict(fit1) - coef(fit1)[1]
# If subset is not NULL, find the subscripts of the cases used:
    if(!is.null(call$subset)){
       rows <- row.names(model.matrix(update(fit1, subset=NULL))) %in% 
               row.names(model.matrix(fit1))
       } else {
       rows <- rep(TRUE, dim(model.matrix(fit1))[1])}
       bign <- length(rows)
       ans <- rep(NA, bign)
       ans[rows] <- p1
       p1 <- ans
# update the formula in cl to fit parallel within group regression
    cl[[2]] <- update(as.formula(formula), as.formula("~ p1*g"))
# update cl$data
    if (is.null(cl$data)){
       cl$data <- data.frame(p1=p1, g=g)} else {
       cl$data <- data.frame(data, p1=p1, g=g)}
# fit parallel lines model  
    fit2 <- eval(cl)
    coef2 <- coef(fit2)
# group.data.subset excludes the subset
    group.data.subset <- 
              model.matrix(fit2$terms, fit2$model, contrasts)[ , 2+1:(l-1)]
    group.data.subset <- data.frame(group=group.data.subset)
# again, fix for non-null subset
    if(!is.null(subset)){
       group.data <- data.frame(matrix(NA, nrow=bign, 
          ncol=dim(group.data.subset)[2]))
       group.data[rows, ] <- group.data.subset} else {
       group.data <- group.data.subset}    
# construct new data with X, y, and G
    formula<-as.formula(formula)
    y.name<-as.character(attr(fit1$terms, "predvars")[2])
    x.name<-attr(fit1$terms, "term.labels")
    g.name<-paste(gname, levels(g)[2:l], sep="")
    names(group.data) <- g.name
    ncols <- length(x.name)+length(y.name)+length(g.name)
    p<-length(x.name)
# generate model
    para.name<-c("eta0", "eta1")
    form1<-paste("eta1 *", x.name[1])
    for(i in 2:p) {
       eta<-paste("eta", i, sep="")
       form1<-paste(form1, "+", eta, "*", x.name[i])
       para.name<-c(para.name, eta)  
    } 
   form2<-paste("eta0 +", form1)
   for(i in 1:(l-1)) {
      G<-g.name[i]
      th0<-paste("th0", i+1, sep="")
      th1<-paste("th1", i+1, sep="")
      form2<-paste(form2, "+", G, "* (", th0, "+", th1, "* (", form1, "))") 
      para.name<-c(para.name, c(th0, th1))
   } 
# New June 6, 2005
   if (is.null(weights))
      {form<-as.formula(paste(y.name, "~", form2))} else
      {wts <- substitute(weights)
       form<-as.formula(paste("sqrt(", wts, ")*", y.name, "~", 
             "sqrt(", wts, ")*(", form2, ")", sep=""))}
# End change
   start <- c(coef2[1], coef2[2]*coef1[-1], coef2[3:(l+1)], 
                coef2[(l+2):(2*l)]/coef2[2])
   names(start) <- c(paste("eta", 0:(length(x.name)), sep=""), 
                     paste("th0", 2:l, sep=""), paste("th1", 2:l, sep=""))
   if (is.null(call$data)){
        for (j in 1:length(g.name)) assign(g.name[j], group.data[ , j])
        obj.nl<-nls(form, start=start, subset=rows, 
            control=control, ...)} else {
        obj.nl<-nls(form, cbind(data, group.data), start=start, 
            subset=rows, control=control, ...)}  
# evaluate linear part, but only for the subset of cases used
   eta.name<-para.name[2:(p+1)]
   eta.value<-coef(obj.nl)[2:(p+1)]
   for(i in 1:p)
      assign(eta.name[i], eta.value[i])
   form1.expr<-parse(text=form1)
   envr <-
     if (is.null(call$data) & is.null(subset)){group.data} else{
      if(is.null(call$data)) group.data.subset else {
       if(is.null(subset)) as.data.frame(cbind(data, group.data)) else{
        as.data.frame(cbind(data[rows, ], group.data.subset))}}}
   if (!is.null(subset)) envr <- as.data.frame(envr[rows, ])
   linear.part<-eval(form1.expr, envr)
   ans <- NULL
   ans$nls.fit <- obj.nl
   ans$linear.part <- linear.part
   ans$group <- if (!is.null(subset)) g[rows] else g
   ans$call <- call
   class(ans)<-c("pod")
   options(op) # warnings ok again
   ans
}}

print.pod<-        function(x, ...){ print(x$nls.fit) }
summary.pod <-     function(object, ...){ summary(object$nls.fit, ...) }
coef.pod <-        function(object, ...){ coef(object$nls.fit, ...)}
deviance.pod <-    function(object, ...){ deviance(object$nls.fit, ...)} 
vcov.pod <-        function(object, ...){ vcov(object$nls.fit, ...)}
residuals.pod <-   function(object, ...){ residuals(object$nls.fit)}
formula.pod<-      function(x, ...){      formula(x$call)}
fitted.pod <-      function(object, ...){ fitted(object$nls.fit, ...)}
df.residual.pod <- function(object, ...){
                                 length(resid(object))-length(coef(object))}
predict.pod     <- function(object, ...){ predict(object$nls.fit)}

# Plot one dimensional models by group
plot.pod.lm <- function(x, colors=1:nlevels(x$group), 
      pch=1:nlevels(x$group), key="topleft", identify=FALSE, 
      xlab="Linear Predictor", ylab=as.character(c(formula(x)[[2]])), ...) {
  mean.function <- x$pod.mean.function
  if(mean.function == "general") stop("No 2D plot for the general pod model")
  g1 <- x$group
  g1.name <- attr(x$group, "name")
  levels(g1) <- 1:nlevels(x$group)
  g1 <- as.numeric(as.character(g1))
  gloc <- match("group", names(x$call))
  gname <- as.character(x$call[[gloc]])
# common regressions
  if(mean.function == "common"){
    plot(predict(x), x$model[ , 1], ylab=ylab, pch=pch[g1], col=colors[g1], 
    xlab=paste(xlab, ", ignore groups", sep=""), ...)
    abline(lm(x$model[ , 1]~predict(x)))}
# parallel regressions
  if(mean.function == "parallel"){
      c2 <- coef(x)
      xp <-0
      for (j in 2:(length(c2)-nlevels(x$group)+1)) 
         xp <- xp +c2[j] * model.matrix(x)[ , j]
      plot(xp, x$model[ , 1], pch=pch[g1], col=colors[g1], 
       ylab=paste(ylab, ", Groups = ", g1.name, sep=""), 
       xlab=paste(xlab,  ", parallel mean function", sep=""), ...)
      for (j in 1:nlevels(x$group))
       abline(if(j==1) c2[1] else c2[1]+c2[length(c2)-nlevels(x$group)+j], 1, 
           lty=j, col=colors[j])}
# key
  if (class(key) == "logical") {
   if (key == TRUE) {
      print("Click mouse on plot to locate the key, or press Escape")
      loc <-locator(n=1) 
      legend(loc[1], loc[2], legend = as.character(levels(x$group)), 
           lty=1:nlevels(x$group), col=colors[1:nlevels(x$group)], 
           pch=pch[1:nlevels(x$group)])}}
   else { 
      loc <- key
      legend(loc[1], loc[2], legend = as.character(levels(x$group)), 
           lty=1:nlevels(x$group), col=colors[1:nlevels(x$group)], 
           inset=0.01, pch=pch[1:nlevels(x$group)])}
# identify
  if(identify == TRUE){
      identify(xp, x$model[ , 1], row.names(x$model))}
  invisible()}

plot.pod <-function(x, colors=1:nlevels(x$group), 
  pch=1:nlevels(x$group), key="topleft", identify=FALSE, 
  xlab="Linear Predictor", ylab=as.character(c(formula(x)[[2]])), ...)
{ 
   yp<-residuals(x)+fitted(x)
   group <- x$group
   gloc <- match("group", names(x$call))
   xp<-x$linear.part
   g1 <- x$group
   levels(g1) <- 1:nlevels(group)
   g1 <- as.numeric(as.character(g1))
   gname <- as.character(x$call[[gloc]])
   plot(xp, yp, pch=pch[g1], col=colors[g1], 
    ylab=paste(ylab, ", Groups = ", as.character(x$call$group), sep=""), 
    xlab=paste(xlab,  ", pod mean function", sep=""), ...)
      for (j in 1:nlevels(group))
        {abline(lm(yp~xp, subset=g1==j), lty=j, col=colors[j])}
   if (class(key) == "logical") {
   if (key == TRUE) {
      print("Click mouse on plot to locate the key, or press Escape")
      loc <-locator(n=1) 
      legend(loc[1], loc[2], legend = as.character(levels(group)), 
           lty=1:nlevels(group), col=colors[1:nlevels(group)], 
           pch=pch[1:nlevels(group)])}} else
    { 
      loc <- key
      legend(loc[1], loc[2], legend = as.character(levels(group)), 
           lty=1:nlevels(group), col=colors[1:nlevels(group)], 
           inset=0.01, pch=pch[1:nlevels(group)])}
   invisible()
}

anova.pod <- 
function (object, scale = 0, test = "F", ...) 
{   
    m1 <- update(object, mean.function="common")
    m2 <- update(object, mean.function="parallel")
    m4 <- update(object, mean.function="general")
    objects <- list(m1, m2, object, m4)
    resdf  <- as.numeric(lapply(objects, df.residual))
    resdev <- as.numeric(lapply(objects, deviance))
    table <- data.frame(resdf, resdev, c(NA, -diff(resdf)), c(NA, 
        -diff(resdev)))
    variables <- c("1: common", "2: parallel", "3: pod", "4: pod + 2fi")
    dimnames(table) <- list(variables, c("Res.Df", "RSS", "Df", 
        "Sum of Sq"))
    title <- paste("POD Analysis of Variance Table for ", 
            deparse(formula(objects[[1]])[[2]]), ", grouped by ", 
            as.character(object$call$group), "\n" , sep="")
    topnote <- c(paste("1: ", deparse(formula(objects[[1]])), sep=""), 
                 paste("2: ", deparse(formula(objects[[2]])), sep=""), 
                 paste("3: ", deparse(formula(object$nls.fit)), sep=""), 
                 paste("4: ", deparse(formula(objects[[4]])), sep=""))
    if (!is.null(test)) {
        bigmodel <- order(resdf)[1]
        scale <- if (scale > 0) 
            scale
        else resdev[bigmodel]/resdf[bigmodel]
        table <- stat.anova(table = table, test = test, scale = scale, 
            df.scale = resdf[bigmodel], n = length(objects[bigmodel$residuals]))
    }
    structure(table, heading = c(title, topnote), class = c("anova", 
        "data.frame"))
}


#######################################################################
# Chapter 7, Transformations
#######################################################################
# Transformations --- all moved to 'car'

####################################################################
#inv.tran.plot now 'invTranPlot' in 'car'

##########################################################################
# replaced by 'powerTransform' in 'car'
    
    
################################################################################
#    Chapter 8
################################################################################
# Test for curvature in a residual plot (Chapter 8)
# Residual plots, and curvature tests in 'car'
  
#resid.curv.test = residCurvTest in 'car'

#resplot = residualPlot in 'car'

#############################################
# marginal model plots    Rev 10/30/07
# mmps = mmps in 'car'

################################################################
# Chapter 9
################################################################

#inf.index = infIndexPlot in 'car'

#outlier.t.test = outlierTest in 'car'
