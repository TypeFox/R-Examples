## nm.mhurdle : extract the names of a relevant subset of the coefficients
## sub.murdle : extract the indexes of a relevant subset of the coefficients
## coef.mhurdle
## vcov.mhurdle
## logLik.mhurdle
## print.mhurdle
## summary.mhurdle
## coef.summary.mhurdle
## print.summary.mhurdle
## fitted.mhurdle
## predict.mhurdle
## update.mhurdle
## rsq

nm.mhurdle <- function(object,
                       which = c("all", "h1", "h2", "h3", "sd", "corr", "tr"),
                       ...){
    coefnames <- object$coef.names
    which <- match.arg(which)
    K <- sapply(coefnames,length)
    if (which == "all"){
        h2.names <- paste("h2", coefnames$h2,sep = ".")
        h1.names <- h3.names <- NULL
        if (! is.null(coefnames$h1)) h1.names <- paste("h1", coefnames$h1,sep = ".")
        if (! is.null(coefnames$h3)) h3.names <- paste("h3", coefnames$h3,sep = ".")
        if (length(coefnames$sd) == 1) sd.names <- "sd"
        else sd.names <- paste("sd", coefnames$sd, sep = ".")
        corr.names <- coefnames$corr
        tr.names <- coefnames$tr
        result <- c(h1.names, h2.names, h3.names, sd.names, corr.names, tr.names)
    }
    else{
        result <- coefnames[[which]]
        if (is.null(result)) stop(paste("no", which, "coefficient\n"))
    }
    result
}

sub.mhurdle <- function(object,
                        which = c("all", "h1", "h2", "h3", "sd", "corr", "tr"),
                        ...){
  # there is no need to check if the coefficient is relevant at it has
  # been checked previously by the nm.mhurdle function
  which <- match.arg(which)
  K <- sapply(object$coef.names, length)
  if (which == "all")  sub <- 1:length(object$coefficients)
  if (which == "h2")   sub <- (K[[1]] + 1):(K[[1]] + K[[2]])
  if (which == "h1")   sub <- 1:K[[1]]
  if (which == "h3")   sub <- (K[[1]] + K[[2]] + 1):(K[[1]] + K[[2]] + K[[3]])
  if (which == "corr") sub <- K[[1]] + K[[2]] + K[[3]] + K[[4]] + 1
  if (which == "sd")   sub <- (K[[1]] + K[[2]] + K[[3]] + 1):(K[[1]] + K[[2]] + K[[3]] + K[[4]])
  if (which == "tr")   sub <- (K[[1]] + K[[2]] + K[[3]] + K[[4]])
  sub
}

coef.mhurdle <- function(object,
                        which = c("all", "h1", "h2", "h3", "sd", "corr", "tr"),
                      ...){
  which <- match.arg(which)
  nm <- nm.mhurdle(object, which)
  sub <- sub.mhurdle(object, which)
  result <- object$coefficients[sub]
  names(result) <- nm
  result
}

vcov.mhurdle <- function(object,
                        which = c("all", "h1", "h2", "h3", "sd", "corr", "tr"),
                      ...){
  which <- match.arg(which)
  nm <- nm.mhurdle(object, which)
  sub <- sub.mhurdle(object, which)
  result <- object$vcov
  result <- result[sub, sub]
  if (is.matrix(result)) rownames(result) <- colnames(result) <- nm
  else names(result) <- nm
  result
}

logLik.mhurdle <- function(object, naive = FALSE, ...){
    if (naive) result <- object$naive$logLik
    else{
        result <- object$logLik
        attrlogLik <- attributes(result)
        result <- sum(object$logLik)
        attributes(result) <- attrlogLik
    }
    result
}


print.mhurdle <- function (x, digits = max(3, getOption("digits") - 2),
                        width = getOption("width"), ...){
  cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
  if (length(coef(x))) {
    cat("Coefficients:\n")
    print.default(format(coef(x), digits = digits), print.gap = 2, 
                  quote = FALSE)
  }
  else cat("No coefficients\n")
  cat("\n")
  invisible(x)
}

summary.mhurdle <- function (object,...){
  b <- coef(object)
  std.err <- sqrt(diag(vcov(object)))
  z <- b/std.err
  p <- 2*(1-pnorm(abs(z)))
  CoefTable <- cbind(b,std.err,z,p)
  colnames(CoefTable) <- c("Estimate","Std. Error","t-value","Pr(>|t|)")
  object$CoefTable <- CoefTable
  object$rsq <- c(coefdet = rsq(object, type = "coefdet"),
                  lratio  = rsq(object, type = "lratio"))
  class(object) <- c("summary.mhurdle","mhurdle")
  return(object)
}


coef.summary.mhurdle <- function(object,
                                 which = c("all", "h1", "h2", "h3", "sd", "corr", "tr"),
                                 ...){
  which <- match.arg(which)
  sub <- sub.mhurdle(object, which)
  nm <- nm.mhurdle(object, which)
  result <- object$CoefTable
  if (!is.null(sub)) result <- result[sub, , drop=FALSE]
  rownames(result) <- nm
  result
}

print.summary.mhurdle <- function(x, digits = max(3, getOption("digits") - 2),
                                  width = getOption("width"), ...){
  cat("\nCall:\n")
  print(x$call)
  cat("\n")
  y <- x$model[,1]
  zeros <- length(y[y==0])/length(y)
  if (zeros>0) cat(paste("Frequency of 0: ",round(zeros,digits=digits),"\n"))
  
  if (!is.null(x$est.stat)){
    cat("\n")
    print(x$est.stat)
  }
  
  cat("\nCoefficients :\n")
  printCoefmat(x$CoefTable,digits=digits)
  cat("\n")
  df <- attr(x$logLik,"df")

  cat(paste("Log-Likelihood: ",
            signif(logLik(x),digits),
            " on ",df," Df\n",sep=""))

  cat("\nR^2 :\n")
  rs <- x$rsq
  cat(paste(" Coefficient of determination :", signif(rs['coefdet'], digits), "\n"))
  cat(paste(" Likelihood ratio index       :", signif(rs['lratio'], digits), "\n"))
  invisible(x)
}

fitted.mhurdle <- function(object, which = c("all", "zero", "positive"), ...){
  which <- match.arg(which)
  switch(which,
         all      = object$fitted.values,
         zero = object$fitted.values[, 1],
         positive = object$fitted.values[, 2]
         )
}

predict.mhurdle <- function(object, newdata = NULL, ...){
    if (is.null(newdata)){
        result <- fitted(object, ...)
    }
    else{
        cl <- object$call
        dist <- ifelse(is.null(cl$dist), TRUE, cl$dist)
        corr <- ifelse(is.null(cl$corr), FALSE, cl$corr)
        m <- model.frame(formula(object), newdata)
        X1 <- model.matrix(formula(object), m, rhs = 1)
        X2 <- model.matrix(formula(object), m, rhs = 2)
        X3 <- model.matrix(formula(object), m, rhs = 3)
        if (length(formula(object))[2] == 4)
            X4 <- model.matrix(formula(object), m, rhs = 4)
        else X4 <- numeric(0)
        y <- model.response(m)
        if (length(X1) == 0) X1 <- NULL
        if (length(X3) == 0) X3 <- NULL
        if (length(X4) == 0) X4 <- NULL
        result <- attr(mhurdle.lnl(coef(object), X1 = X1, X2 = X2, X3 = X3, X4 = X4, y = y,
                                   gradient = FALSE, fitted = TRUE,
                                   dist = dist, corr = corr), "fitted")
    }
    result
}

## a simple copy from mlogit. update with formula doesn't work
## otherwise ????
update.mhurdle <- function (object, new, ...){
  call <- object$call
  if (is.null(call))
    stop("need an object with call component")
  extras <- match.call(expand.dots = FALSE)$...
  if (!missing(new))
    call$formula <- update(formula(object), new)
  if(length(extras) > 0) {
    existing <- !is.na(match(names(extras), names(call)))
    ## do these individually to allow NULL to remove entries.
    for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
    if(any(!existing)) {
      call <- c(as.list(call), extras[!existing])
      call <- as.call(call)
    }
  }
  for (i in 1:length(attr(call$formula, "rhs"))){
    # update.Formula returns "1 - 1" instead of 0 for empty parts
    zero <- paste(deparse(attr(call$formula, "rhs")[[i]])) == as.character("1 - 1")
    if (zero) attr(call$formula, "rhs")[[i]] <- 0
  }
  eval(call, parent.frame())
}

rsq <- function(object,
                type = c("coefdet", "lratio"),
                adj = FALSE,
                r2pos=c("rss","ess","cor")){
    
    type <- match.arg(type)
    r2pos <- match.arg(r2pos)
    
    K1 <- length(object$coef.names$h1)
    K2 <- length(object$coef.names$h2)
    K3 <- length(object$coef.names$h3)
    
    y <- model.response(model.frame(object))
    n <- length(y)
    no <- sum(y == 0)
    po <- mean(y == 0)
    pp <- 1 - po

    K <- length(coef(object))
    Ko <- length(object$naive$coefficients)
    
    if (type == "lratio"){
        if (!adj) R2 <- 1 - logLik(object) / logLik(object, naive = TRUE)
        else R2 <- 1 - (logLik(object) - K) / (logLik(object, naive = TRUE) - Ko)
        R2 <- as.numeric(R2)
    }
    if (type == "coefdet"){
        ym <- mean(y)
        yf <- fitted(object, "positive") * (1 - fitted(object, "zero"))
        R2 <- switch(r2pos,
                     ess = ifelse(adj,
                         sum( (yf - ym) ^ 2) / sum( (y - ym) ^ 2) * (n - K) / (n - Ko),
                         sum( (yf - ym) ^ 2) / sum( (y - ym) ^ 2)),
                     rss = ifelse(adj,
                         1 - (n - Ko) / (n - K) * sum( (y - yf) ^ 2) / sum( (y - ym) ^ 2),
                         1 - sum( (y - yf) ^ 2) / sum( (y - ym) ^ 2)),
                     cor = ifelse(adj,
                         stop("no adjusted R2 using the correlation formula"),
                         cor(y, yf) ^ 2
                         )
                     )
    }
    R2
}
