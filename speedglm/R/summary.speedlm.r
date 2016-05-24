summary.speedlm <- function (object, correlation = FALSE, ...) 
{
  if (!inherits(object, "speedlm")) 
    stop("object must be an object of class speedlm")
  z <- object
  n <- if (is.null(z$weights)) z$nobs else z$nobs - z$zero.w
  nvar <- z$nvar
  rdf <- z$df.residual
  if (z$method=='qr') {
    z$XTX <- z$XTX[z$ok,z$ok]
  }
  var_res <- as.numeric(z$RSS)/rdf
  se_coef <- rep(NA, z$nvar)
  inv <- solve.qr(qr(z$XTX))
  se_coef[z$ok] <- sqrt(var_res * diag(inv))
  t1 <- z$coefficients/se_coef
  p <- 2 * pt(abs(t1), df = z$df.residual, lower.tail = FALSE)
  ip <- !is.na(p)
  if (is.null(z$weights)) {
    X1X <- z$X1X[z$ok]
    if (z$intercept) {
      X1X <- matrix(kronecker(X1X, X1X), z$rank, z$rank)
      mss <- crossprod(z$coef, z$XTX - X1X/n) %*% z$coef
    } else mss <- crossprod(z$coef, z$XTX) %*% z$coef
  } else {
    XW1 <- z$XW1[z$ok]
    mss <- if (z$intercept) {
      XWX <- matrix(kronecker(XW1, XW1), length(XW1), length(XW1))
      XW1 <- matrix(kronecker(XW1 * z$SW, XW1), z$rank, z$rank, byrow = TRUE)
      crossprod(z$coef, z$XTX - 2 * XWX/z$SW + XW1/(z$SW^2)) %*% z$coef
    } else crossprod(z$coef, z$XTX) %*% z$coef
  }
  rss <- z$RSS
  if (nvar != (z$intercept)) {
    df.int <- if (z$intercept) 1L else 0L
    r.squared <- as.numeric(mss/(mss + rss))
    adj.r.squared <- 1 - (1 - r.squared) * ((n - df.int)/rdf)
    fstatistic <- c(value = (as(mss, "numeric")/(z$rank - 
                                                   df.int))/var_res, numdf = z$rank - df.int, dendf = rdf)
    f.pvalue <- 1 - pf(fstatistic[1], fstatistic[2], fstatistic[3])
  }
  else {
    fstatistic <- f.pvalue <- NULL
    r.squared <- adj.r.squared <- 0
  }  
  param <- data.frame(coef = z$coefficients, se = se_coef, 
                      t = t1, `p-value` = p)
  keep <- match(c("call", "terms", "frame", "ok", "RSS", "rank"), 
                names(object), 0)
  ans <- c(object[keep], list(coefficients = param, var.res = var_res, 
                              df.residuals = z$df.residual, nobs = z$nobs, r.squared = r.squared, 
                              adj.r.squared = adj.r.squared, correlation = correlation, 
                              fstatistic = fstatistic, f.pvalue = f.pvalue, rdf = rdf, 
                              cov.scaled = inv, intercept = (nvar != (z$intercept))))
  if (correlation) 
    ans$correl <- inv * as.numeric(var_res)/outer(se_coef[z$ok], 
                                                  se_coef[z$ok])
  class(ans) <- "summary.speedlm"
  return(ans)
}

print.summary.speedlm <- function(x,digits=max(3,getOption("digits")-3),...){
  
  x$coefficients$coef <- if (any(abs(na.omit(x$coefficients$coef))<0.0001))
    format(x$coefficients$coef,scientific=TRUE,
           digits=4) else  round(x$coefficients$coef,digits=6)
  x$coefficients$se <- if (any(na.omit(x$coefficients$se)<0.0001))
    format(x$coefficients$se,scientific=TRUE,digits=4) else
      round(x$coefficients$se,digits=6)
  x$coefficients$t <- round(x$coefficients$t,digits=3)
  x$coefficients$p.value <- if (any(na.omit(x$coefficients$p.value)<0.0001))
    format(x$coefficients$p.value,scientific=TRUE,
           digits=4) else round(x$coefficients$p.value,
                                digits=6)
  s<-sum(Vectorize(is.na(x$coefficients$coef)))
  cat("Linear Regression Model of class 'speedlm':\n")  
  if (!is.null(x$call)) cat("\nCall: ", deparse(x$call), "\n\n")  
  if (length(x$coef)) {
    cat("Coefficients:\n")
    cat(" ------------------------------------------------------------------", "\n")
    sig <- function(z) {
      if (z < 0.001) 
        "***"
      else if (z < 0.01) 
        "** "
      else if (z < 0.05) 
        "*  "
      else if (z < 0.1) 
        ".  "
      else "   "
    }
    sig.1 <- sapply(x$coefficients$p.value, sig)
    est.1 <- cbind(format(x$coefficients, digits = digits), sig.1)
    colnames(est.1)[ncol(est.1)] <- ""
    print(est.1)
    cat("\n")
    cat("-------------------------------------------------------------------", 
        "\n")
    cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1", 
        "\n") 
  } else cat("No coefficients\n")
  cat("---\n")
  if (x$intercept)
    cat("Residual standard error: ",round(sqrt(x$var.res),6)," on ",x$rdf,
        " degrees of freedom;\n",
        "observations: ",x$nobs,";  R^2: ",format(x$r.squared,digits=3),
        ";  adjusted R^2: ",format(x$adj.r.squared,digits=3),";\n",
        "F-statistic: ",format(x$fstatistic[1],digits=4)," on ", x$fstatistic[2],
        " and ",x$fstatistic[3]," df;  p-value: ",format(x$f.pvalue,digits=6),
        ".\n",sep="") else 
          cat("Residual standard error: ",round(sqrt(x$var.res),6)," on ",x$rdf,
              " degrees of freedom\n",sep="")   
  if (s==1) cat("One coefficient not defined because of singularities. \n")
  if (s>1) cat(s," coefficients not defined because of singularities. \n")
  invisible(x)
  if (x$correlation) {
    cat("---\n")
    cat("Correlation of Coefficients:\n")
    x$correl[upper.tri(x$correl,diag=TRUE)]<-NA
    print(x$correl[-1,-nrow(x$correl)],na.print="",digits=2)
  }
}



coef.speedlm <- function(object,...)  object$coefficients


vcov.speedlm <- function(object,...){
  z <- object
  var_res <- z$RSS/z$df.residual
  rval<- var_res * solve(z$XTX)
  rval
}

logLik.speedlm <- function(object,...){
    p <- object$rank
    N <- object$nobs
    if (is.null(object$pw)) pw <- 1 else  {
      N <- object$nobs - object$zero.w
      pw <- object$pw
    }  
    N0 <- N
    val <- 0.5 * (log(pw) - N * (log(2 * pi) + 1 - log(N) +  log(object$RSS)))
    attr(val, "nall") <- N0
    attr(val, "nobs") <- N
    attr(val, "df") <- p + 1
    class(val)<-"logLik.speedlm"
    val
}

print.logLik.speedlm <- function(x, digits = getOption("digits"), ...) {
    cat("'log Lik.' ", paste(format(c(x), digits = digits), collapse = ", "), 
        " (df=", format(attr(x, "df")), ")\n", sep = "")
    invisible(x)
}

AIC.speedlm<-function (object,...,k = 2){ 
    p <- object$rank
    N <- object$nobs
    if (is.null(object$pw)) pw <- 1 else  {
      N <- object$nobs - object$zero.w
      pw <- object$pw
    }  
    val <- -(log(pw) - N * (log(2 * pi) + 1 - log(N) + log(object$RSS)))+k*(p+1)
    val
}



print.speedlm <- function(x,digits = max(3, getOption("digits") - 3),...){
  cat("Linear Regression Model of class 'speedlm':\n")
  if (!is.null(x$call)) cat("\nCall: ", deparse(x$call), "\n\n")
  if (length(x$coef)) {
        cat("Coefficients:\n")
        print.default(format(x$coefficients, digits = digits), print.gap = 2,
            quote = FALSE)
    } else cat("No coefficients\n")
    cat("\n")
    invisible(x)
}


extractAIC.speedlm<-function(fit, scale=0, k=2,...) 
{
  n <- fit$nobs
  edf <- n - fit$df.residual
  RSS <- fit$RSS
  dev <- if (scale > 0) 
    RSS/scale - n
  else n * log(RSS/n)
  c(edf, dev + k * edf)
}

drop1.speedlm<-function (object, scope, scale = 0, all.cols = TRUE,test = c("none","Chisq", "F"), k = 2, data, ...) {
  #	check_exact(object)
  x <- model.matrix(object$terms, (if(missing(data)) object$model else data))  
  offset <- model.offset(model.frame(object))
  iswt <- !is.null(wt <- object$weights)
  n <- nrow(x)
  asgn <- attr(x, "assign")
  tl <- attr(object$terms, "term.labels")
  if (missing(scope)) 
    scope <- drop.scope(object)
  else {
    if (!is.character(scope)) 
      scope <- attr(terms(update.formula(object, scope)), 
                    "term.labels")
    if (!all(match(scope, tl, 0L) > 0L)) 
      stop("scope is not a subset of term labels")
  }
  ndrop <- match(scope, tl)
  ns <- length(scope)
  rdf <- object$df.residual
  chisq <- object$RSS #deviance.lm(object)
  dfs <- numeric(ns)
  RSS <- numeric(ns)
  y <- (if(missing(data)) object$model else data)[,rownames(attr(object$terms,'factors'))[attr(object$terms,'response')]]#object$residuals + object$fitted.values
  na.coef <- seq_along(object$coefficients)[!is.na(object$coefficients)]
  for (i in seq_len(ns)) {
    ii <- seq_along(asgn)[asgn == ndrop[i]]
    jj <- setdiff(if (all.cols) 
      seq(ncol(x))
      else na.coef, ii)
    z <- if (iswt) 
      speedlm.wfit(y,x[, jj, drop = FALSE], wt, offset = offset,...)
    else speedlm.fit(y,x[, jj, drop = FALSE], offset = offset,method=object$method,...)
    dfs[i] <- z$rank
    RSS[i] <- z$RSS
  }
  scope <- c("<none>", scope)
  dfs <- c(object$rank, dfs)
  RSS <- c(chisq, RSS)
  if (scale > 0) 
    aic <- RSS/scale - n + k * dfs
  else aic <- n * log(RSS/n) + k * dfs
  dfs <- dfs[1L] - dfs
  dfs[1L] <- NA
  aod <- data.frame(Df = dfs, `Sum of Sq` = c(NA, RSS[-1L] - 
                                                RSS[1L]), RSS = RSS, AIC = aic, row.names = scope, check.names = FALSE)
  if (scale > 0) 
    names(aod) <- c("Df", "Sum of Sq", "RSS", "Cp")
  test <- match.arg(test)
  if (test == "Chisq") {
    dev <- aod$"Sum of Sq"
    if (scale == 0) {
      dev <- n * log(RSS/n)
      dev <- dev - dev[1L]
      dev[1L] <- NA
    }
    else dev <- dev/scale
    df <- aod$Df
    nas <- !is.na(df)
    dev[nas] <- safe_pchisq(dev[nas], df[nas], lower.tail = FALSE)
    aod[, "Pr(>Chi)"] <- dev
  }
  else if (test == "F") {
    dev <- aod$"Sum of Sq"
    dfs <- aod$Df
    rdf <- object$df.residual
    rms <- aod$RSS[1L]/rdf
    Fs <- (dev/dfs)/rms
    Fs[dfs < 1e-04] <- NA
    P <- Fs
    nas <- !is.na(Fs)
    P[nas] <- safe_pf(Fs[nas], dfs[nas], rdf, lower.tail = FALSE)
    aod[, c("F value", "Pr(>F)")] <- list(Fs, P)
  }
  head <- c("Single term deletions", "\nModel:", deparse(formula(object)), 
            if (scale > 0) paste("\nscale: ", format(scale), "\n"))
  class(aod) <- c("anova", "data.frame")
  attr(aod, "heading") <- head
  aod
}


add1.speedlm<-function (object, scope, scale = 0, test = c("none", "Chisq","F"), x = NULL, k = 2,data, ...) 
{
  Fstat <- function(table, RSS, rdf) {
    dev <- table$"Sum of Sq"
    df <- table$Df
    rms <- (RSS - dev)/(rdf - df)
    Fs <- (dev/df)/rms
    Fs[df < .Machine$double.eps] <- NA
    P <- Fs
    nnas <- !is.na(Fs)
    P[nnas] <- safe_pf(Fs[nnas], df[nnas], rdf - df[nnas], 
                       lower.tail = FALSE)
    list(Fs = Fs, P = P)
  }
  #	check_exact(object)
  if (missing(scope) || is.null(scope)) 
    stop("no terms in scope")
  if (!is.character(scope)) 
    scope <- add.scope(object, update.formula(object, scope))
  if (!length(scope)) 
    stop("no terms in scope for adding to object")
  oTerms <- attr(object$terms, "term.labels")
  int <- attr(object$terms, "intercept")
  ns <- length(scope)
  y <- (if(missing(data)) object$model else data)[,rownames(attr(object$terms,'factors'))[attr(object$terms,'response')]]#object$residuals + object$fitted.values
  dfs <- numeric(ns + 1)
  RSS <- numeric(ns + 1)
  names(dfs) <- names(RSS) <- c("<none>", scope)
  add.rhs <- paste(scope, collapse = "+")
  add.rhs <- eval(parse(text = paste("~ . +", add.rhs)))
  new.form <- update.formula(object, add.rhs)
  Terms <- terms(new.form)
  if (is.null(x)) {
    fc <- object$call
    fc$formula <- Terms
    fob <- list(call = fc, terms = Terms)
    class(fob) <- 'speedlm' # oldClass(object)
    m <- model.frame(fob, (if(missing(data)) object$model else data))
    x <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
    offset <- model.offset(m)
    wt <- model.weights(m)
    oldn <- length(y)
    y <- model.response(m, "numeric")
    newn <- length(y)
    if (newn < oldn) 
      warning(sprintf(ngettext(newn, "using the %d/%d row from a combined fit", 
                               "using the %d/%d rows from a combined fit"), 
                      newn, oldn), domain = NA)
  }
  else {
    wt <- object$weights
    offset <- object$offset
  }
  n <- nrow(x)
  Terms <- attr(Terms, "term.labels")
  asgn <- attr(x, "assign")
  ousex <- match(asgn, match(oTerms, Terms), 0L) > 0L
  if (int) 
    ousex[1L] <- TRUE
  iswt <- !is.null(wt)
  X <- x[, ousex, drop = FALSE]
  z <- if (iswt) 
    speedlm.wfit(y,X, wt, offset = offset,...)
  else speedlm.fit(y, X, offset = offset, method=object$method,...)
  dfs[1L] <- z$rank
  class(z) <- "speedlm"
  RSS[1L] <- z$RSS #deviance(z)
  sTerms <- sapply(strsplit(Terms, ":", fixed = TRUE), function(x) paste(sort(x), 
                                                                         collapse = ":"))
  for (tt in scope) {
    stt <- paste(sort(strsplit(tt, ":")[[1L]]), collapse = ":")
    usex <- match(asgn, match(stt, sTerms), 0L) > 0L
    X <- x[, usex | ousex, drop = FALSE]
    z <- if (iswt) 
      speedlm.wfit(y, X, wt, offset = offset,...)
    else speedlm.fit(y, X, offset = offset,method=object$method,...)
    class(z) <- "speedlm"
    dfs[tt] <- z$rank
    RSS[tt] <- z$RSS #deviance(z)
  }
  if (scale > 0) 
    aic <- RSS/scale - n + k * dfs
  else aic <- n * log(RSS/n) + k * dfs
  dfs <- dfs - dfs[1L]
  dfs[1L] <- NA
  aod <- data.frame(Df = dfs, `Sum of Sq` = c(NA, RSS[1L] - 
                                                RSS[-1L]), RSS = RSS, AIC = aic, row.names = names(dfs), 
                    check.names = FALSE)
  if (scale > 0) 
    names(aod) <- c("Df", "Sum of Sq", "RSS", "Cp")
  test <- match.arg(test)
  if (test == "Chisq") {
    dev <- aod$"Sum of Sq"
    if (scale == 0) {
      dev <- n * log(RSS/n)
      dev <- dev[1L] - dev
      dev[1L] <- NA
    }
    else dev <- dev/scale
    df <- aod$Df
    nas <- !is.na(df)
    dev[nas] <- safe_pchisq(dev[nas], df[nas], lower.tail = FALSE)
    aod[, "Pr(>Chi)"] <- dev
  }
  else if (test == "F") {
    rdf <- object$df.residual
    aod[, c("F value", "Pr(>F)")] <- Fstat(aod, aod$RSS[1L], 
                                           rdf)
  }
  head <- c("Single term additions", "\nModel:", deparse(formula(object)), 
            if (scale > 0) paste("\nscale: ", format(scale), "\n"))
  class(aod) <- c("anova", "data.frame")
  attr(aod, "heading") <- head
  aod
}



nobs.speedlm <- function(object, use.fallback = FALSE,...) 
  if (!is.null(w <- object$weights)) sum(w != 0) else object$nobs




