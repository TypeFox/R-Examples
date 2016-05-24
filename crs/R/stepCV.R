# file MASS/R/stepAIC.R
# copyright (C) 1994-2007 W. N. Venables and B. D. Ripley
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#

## The file add1.lm.cv is a modified version of add1.lm found in
## add.R. To the original object I append the cross-validation
## score. I do this so that I can modify stepAIC to admit the
## cross-validation score rather than AIC/BIC (k=log(n)).

extractCV <- function(model,...) {
  ## Maximizing hence negative sign
  c(model$rank, mean(residuals(model)^2/(1-hatvalues(model))^2))
}


check_exact <- function(object)
{
    w <- object$weights
    if(is.null(w)) {
        mss <- sum(object$fitted.values^2)
        rss <- sum(object$residuals^2)
    } else {
        mss <- sum(w * object$fitted.values^2)
        rss <- sum(w * object$residuals^2)
    }
    if(rss < 1e-10*mss)
        warning("attempting model selection on an essentially perfect fit is nonsense", call. = FALSE)
}
add1.lm.cv <- function (object, scope, scale = 0, test = c("none", "Chisq", 
    "F"), x = NULL, k = 2, ...) 
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
    check_exact(object)
    if (missing(scope) || is.null(scope)) 
        stop("no terms in scope")
    if (!is.character(scope)) 
        scope <- add.scope(object, update.formula(object, scope))
    if (!length(scope)) 
        stop("no terms in scope for adding to object")
    oTerms <- attr(object$terms, "term.labels")
    int <- attr(object$terms, "intercept")
    ns <- length(scope)
    y <- object$residuals + object$fitted.values
    dfs <- numeric(ns + 1)
    RSS <- numeric(ns + 1)
    CV <- numeric(ns + 1) # added jracine
    names(CV) <- names(dfs) <- names(RSS) <- c("<none>", scope) # added term 1 jracine
    add.rhs <- paste(scope, collapse = "+")
    add.rhs <- eval(parse(text = paste("~ . +", add.rhs)))
    new.form <- update.formula(object, add.rhs)
    Terms <- terms(new.form)
    if (is.null(x)) {
        fc <- object$call
        fc$formula <- Terms
        fob <- list(call = fc, terms = Terms)
        class(fob) <- oldClass(object)
        m <- model.frame(fob, xlev = object$xlevels)
        x <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
        offset <- model.offset(m)
        wt <- model.weights(m)
        oldn <- length(y)
        y <- model.response(m, "numeric")
        newn <- length(y)
        if (newn < oldn) 
            warning(gettextf("using the %d/%d rows from a combined fit", 
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
      lm.wfit(X, y, wt, offset = offset)
    else lm.fit(X, y, offset = offset)
    dfs[1L] <- z$rank
    class(z) <- "lm"
    RSS[1L] <- deviance(z)
    model <- lm(y~X) # added jracine
    CV[1L] <- mean(residuals(model)^2/(1-hatvalues(model))^2) # added jracine
    sTerms <- sapply(strsplit(Terms, ":", fixed = TRUE), function(x) paste(sort(x), 
        collapse = ":"))
    for (tt in scope) {
        stt <- paste(sort(strsplit(tt, ":")[[1L]]), collapse = ":")
        usex <- match(asgn, match(stt, sTerms), 0L) > 0L
        X <- x[, usex | ousex, drop = FALSE]
        z <- if (iswt) 
            lm.wfit(X, y, wt, offset = offset)
        else lm.fit(X, y, offset = offset)
        class(z) <- "lm"
        dfs[tt] <- z$rank
        RSS[tt] <- deviance(z)
        model <- lm(y~X) # added jracine
        CV[tt] <- mean(residuals(model)^2/(1-hatvalues(model))^2) # added jracine
    }
    if (scale > 0) 
        aic <- RSS/scale - n + k * dfs
    else aic <- n * log(RSS/n) + k * dfs
    dfs <- dfs - dfs[1L]
    dfs[1L] <- NA
    ## Can add here
    aod <- data.frame(Df = dfs, `Sum of Sq` = c(NA, RSS[1L] - 
        RSS[-1L]), RSS = RSS, AIC = aic, CV=CV, row.names = names(dfs), 
        check.names = FALSE) # added CV term jracine
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
        aod[, "Pr(Chi)"] <- dev
    }
    else if (test == "F") {
        rdf <- object$df.residual
        aod[, c("F value", "Pr(F)")] <- Fstat(aod, aod$RSS[1L], 
            rdf)
    }
    head <- c("Single term additions", "\nModel:", deparse(as.vector(formula(object))), 
        if (scale > 0) paste("\nscale: ", format(scale), "\n"))
    class(aod) <- c("anova", "data.frame")
    attr(aod, "heading") <- head
    aod
}

stepCV <-
  function(object, scope, scale = 0,
           direction = c("both", "backward", "forward"),
           trace = 1, keep = NULL, steps = 1000, use.start = FALSE, k = 2, ...)
{
    mydeviance <- function(x, ...)
    {
        dev <- deviance(x)
        if(!is.null(dev)) dev else extractCV(x, k=0)[2L]
    }

    cut.string <- function(string)
    {
        if(length(string) > 1L)
            string[-1L] <- paste("\n", string[-1L], sep = "")
        string
    }

    re.arrange <- function(keep)
    {
        namr <- names(k1 <- keep[[1L]])
        namc <- names(keep)
        nc <- length(keep)
        nr <- length(k1)
        array(unlist(keep, recursive = FALSE), c(nr, nc), list(namr, namc))
    }

    step.results <- function(models, fit, object, usingCp=FALSE)
    {
        change <- sapply(models, "[[", "change")
        rd <- sapply(models, "[[", "deviance")
        dd <- c(NA, abs(diff(rd)))
        rdf <- sapply(models, "[[", "df.resid")
        ddf <- c(NA, abs(diff(rdf)))
        CV <- sapply(models, "[[", "CV")
        heading <- c("Stepwise Model Path \nAnalysis of Deviance Table",
                     "\nInitial Model:", deparse(as.vector(formula(object))),
                     "\nFinal Model:", deparse(as.vector(formula(fit))),
                     "\n")
        aod <-
            if(usingCp)
                data.frame(Step = change, Df = ddf, Deviance = dd,
                           "Resid. Df" = rdf, "Resid. Dev" = rd,
                           Cp = CV, check.names = FALSE)
            else data.frame(Step = change, Df = ddf, Deviance = dd,
                            "Resid. Df" = rdf, "Resid. Dev" = rd,
                            CV = CV, check.names = FALSE)
        attr(aod, "heading") <- heading
        class(aod) <- c("Anova", "data.frame")
        fit$anova <- aod
        fit
    }

    Terms <- terms(object)
    object$formula <- Terms
    if(inherits(object, "lme")) object$call$fixed <- Terms
    else if(inherits(object, "gls")) object$call$model <- Terms
    else object$call$formula <- Terms
    if(use.start) warning("'use.start' cannot be used with R's version of glm")
    md <- missing(direction)
    direction <- match.arg(direction)
    backward <- direction == "both" | direction == "backward"
    forward <- direction == "both" | direction == "forward"
    if(missing(scope)) {
        fdrop <- numeric(0)
        fadd <- attr(Terms, "factors")
        if(md) forward <- FALSE
    } else {
        if(is.list(scope)) {
            fdrop <- if(!is.null(fdrop <- scope$lower))
                attr(terms(update.formula(object, fdrop)), "factors")
            else numeric(0)
            fadd <- if(!is.null(fadd <- scope$upper))
                attr(terms(update.formula(object, fadd)), "factors")
        } else {
            fadd <- if(!is.null(fadd <- scope))
                attr(terms(update.formula(object, scope)), "factors")
            fdrop <- numeric(0)
        }
    }
    models <- vector("list", steps)
    if(!is.null(keep)) keep.list <- vector("list", steps)
    ## watch out for partial matching here.
    if(is.list(object) && (nmm <- match("nobs", names(object), 0)) > 0)
        n <- object[[nmm]]
    else n <- length(residuals(object))
    fit <- object
    bCV <- extractCV(fit, scale, k = k, ...)
    edf <- bCV[1L]
    bCV <- bCV[2L]
    if(is.na(bCV))
        stop("CV is not defined for this model, so stepCV cannot proceed")
    nm <- 1
    Terms <- terms(fit)
    if(trace) {
        cat("Start:  CV=", format(bCV), "\n",
            cut.string(deparse(as.vector(formula(fit)))), "\n\n", sep='')
  utils::flush.console()
    }
    models[[nm]] <- list(deviance = mydeviance(fit), df.resid = n - edf,
                         change = "", CV = bCV)
    if(!is.null(keep)) keep.list[[nm]] <- keep(fit, bCV)
    usingCp <- FALSE
    while(steps > 0) {
        steps <- steps - 1
        CV <- bCV
        ffac <- attr(Terms, "factors")
        ## don't drop strata terms
        if(!is.null(sp <- attr(Terms, "specials")) &&
           !is.null(st <- sp$strata)) ffac <- ffac[-st,]
        scope <- factor.scope(ffac, list(add = fadd, drop = fdrop))
        aod <- NULL
        change <- NULL
        if(backward && length(scope$drop)) {
            aod <- dropterm(fit, scope$drop, scale = scale,
                            trace = max(0, trace - 1), k = k, ...)
            rn <- row.names(aod)
            row.names(aod) <- c(rn[1L], paste("-", rn[-1L], sep=" "))
            ## drop all zero df terms first.
            if(any(aod$Df == 0, na.rm=TRUE)) {
                zdf <- aod$Df == 0 & !is.na(aod$Df)
                nc <- match(c("Cp", "CV"), names(aod))
                nc <- nc[!is.na(nc)][1L]
                ch <- abs(aod[zdf, nc] - aod[1, nc]) > 0.01
                if(any(ch)) {
                    warning("0 df terms are changing CV")
                    zdf <- zdf[!ch]
                }
                ## drop zero df terms first: one at time since they
                ## may mask each other
                if(length(zdf) > 0L)
                    change <- rev(rownames(aod)[zdf])[1L]
            }
        }
        if(is.null(change)) {
            if(forward && length(scope$add)) {
                aodf <- addterm(fit, scope$add, scale = scale,
                                trace = max(0, trace - 1), k = k, ...)
                rn <- row.names(aodf)
                row.names(aodf) <- c(rn[1L], paste("+", rn[-1L], sep=" "))
                aod <-
                    if(is.null(aod)) aodf
                    else rbind(aod, aodf[-1, , drop=FALSE])
            }
            attr(aod, "heading") <- NULL
            if(is.null(aod) || ncol(aod) == 0) break
            ## need to remove any terms with zero df from consideration
            nzdf <- if(!is.null(aod$Df)) aod$Df != 0 | is.na(aod$Df)
            aod <- aod[nzdf, ]
            if(is.null(aod) || ncol(aod) == 0) break
            nc <- match(c("Cp", "CV"), names(aod))
            nc <- nc[!is.na(nc)][1L]
            o <- order(aod[, nc])
            if(trace) {
    print(aod[o,  ])
    utils::flush.console()
      }
            if(o[1L] == 1) break
            change <- rownames(aod)[o[1L]]
        }
        usingCp <- match("Cp", names(aod), 0) > 0
        ## may need to look for a 'data' argument in parent
        fit <- update(fit, paste("~ .", change), evaluate = FALSE)
        fit <- eval.parent(fit)
        if(is.list(fit) && (nmm <- match("nobs", names(fit), 0)) > 0)
            nnew <- fit[[nmm]]
        else nnew <- length(residuals(fit))
        if(nnew != n)
            stop("number of rows in use has changed: remove missing values?")
        Terms <- terms(fit)
        bCV <- extractCV(fit, scale, k = k, ...)
        edf <- bCV[1L]
        bCV <- bCV[2L]
        if(trace) {
            cat("\nStep:  CV=", format(bCV), "\n",
                cut.string(deparse(as.vector(formula(fit)))), "\n\n", sep='')
      utils::flush.console()
  }
        ## add a tolerance as dropping 0-df terms might increase CV slightly
        if(bCV >= CV + 1e-7) break
        nm <- nm + 1
        models[[nm]] <-
            list(deviance = mydeviance(fit), df.resid = n - edf,
                 change = change, CV = bCV)
        if(!is.null(keep)) keep.list[[nm]] <- keep(fit, bCV)
    }
    if(!is.null(keep)) fit$keep <- re.arrange(keep.list[seq(nm)])
    step.results(models = models[seq(nm)], fit, object, usingCp)
}

extractCV.loglm <- function(fit, scale, k = 2, ...)
{
    edf <- fit$n - fit$df
    c(edf,  fit$deviance + k * edf)
}

extractCV.lme <- function(fit, scale, k = 2, ...)
{
    if(fit$method != "ML") stop("CV undefined for REML fit")
    res <- logLik(fit)
    edf <- attr(res, "df")
    c(edf,  -2*res + k * edf)
}

extractCV.gls <- function(fit, scale, k = 2, ...)
{
    if(fit$method != "ML") stop("CV undefined for REML fit")
    res <- logLik(fit)
    edf <- attr(res, "df")
    c(edf,  -2*res + k * edf)
}

terms.gls <- terms.lme <- function(x, ...) terms(formula(x), ...)

# file MASS/R/add.R
# copyright (C) 1994-2008 W. N. Venables and B. D. Ripley
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#

## version to return NA for df = 0, as R did before 2.7.0
safe_pchisq <- function(q, df, ...)
{
    df[df <= 0] <- NA
    pchisq(q=q, df=df, ...)
}
## and to avoid a warning
safe_pf <- function(q, df1, ...)
{
    df1[df1 <= 0] <- NA
    pf(q=q, df1=df1, ...)
}


addterm <-
    function(object, ...) UseMethod("addterm")

addterm.default <-
    function(object, scope, scale = 0, test = c("none", "Chisq"),
             k = 2, sorted = FALSE, trace = FALSE, ...)
{
    if(missing(scope) || is.null(scope)) stop("no terms in scope")
    if(!is.character(scope))
        scope <- add.scope(object, update.formula(object, scope))
    if(!length(scope))
        stop("no terms in scope for adding to object")
#     newform <- update.formula(object,
#                               paste(". ~ . +", paste(scope, collapse="+")))
#     data <- model.frame(update(object, newform)) # remove NAs
#     object <- update(object, data = data)
    ns <- length(scope)
    ans <- matrix(nrow = ns + 1L, ncol = 2L,
                  dimnames = list(c("<none>", scope), c("df", "CV")))
    ans[1L,  ] <- extractCV(object, scale, k = k, ...)
    n0 <- length(object$residuals)
    env <- environment(formula(object))
    for(i in seq(ns)) {
        tt <- scope[i]
        if(trace) {
      message("trying +", tt)
      utils::flush.console()
        }
        nfit <- update(object, as.formula(paste("~ . +", tt)),
                       evaluate = FALSE)
  nfit <- eval(nfit, envir=env) # was  eval.parent(nfit)
  ans[i+1L, ] <- extractCV(nfit, scale, k = k, ...)
        if(length(nfit$residuals) != n0)
            stop("number of rows in use has changed: remove missing values?")
    }
    dfs <- ans[, 1L] - ans[1L, 1L]
    dfs[1L] <- NA
    aod <- data.frame(Df = dfs, CV = ans[, 2L])
    o <- if(sorted) order(aod$CV) else seq_along(aod$CV)
    test <- match.arg(test)
    if(test == "Chisq") {
  dev <- ans[, 2L] - k*ans[, 1L]
  dev <- dev[1L] - dev; dev[1L] <- NA
  nas <- !is.na(dev)
  P <- dev
  P[nas] <- safe_pchisq(dev[nas], dfs[nas], lower.tail=FALSE)
  aod[, c("LRT", "Pr(Chi)")] <- list(dev, P)
    }
    aod <- aod[o, ]
    head <- c("Single term additions", "\nModel:",
              deparse(as.vector(formula(object))))
    if(scale > 0)
        head <- c(head, paste("\nscale: ", format(scale), "\n"))
    class(aod) <- c("anova", "data.frame")
    attr(aod, "heading") <- head
    aod
}

addterm.lm <-
  function(object, scope, scale = 0, test = c("none", "Chisq", "F"),
           k = 2, sorted = FALSE, ...)
{
    Fstat <- function(table, RSS, rdf) {
        dev <- table$"Sum of Sq"
        df <- table$Df
        rms <- (RSS - dev)/(rdf - df)
        Fs <- (dev/df)/rms
        Fs[df < 1e-4] <- NA
        P <- Fs
        nnas <- !is.na(Fs)
  P[nnas] <- pf(Fs[nnas], df[nnas], rdf - df[nnas], lower.tail=FALSE)
        list(Fs=Fs, P=P)
    }

    if(missing(scope) || is.null(scope)) stop("no terms in scope")
#    aod <- stats:::add1.lm(object, scope=scope, scale=scale)[ , -4L]
    aod <- add1.lm.cv(object, scope=scope, scale=scale)[ , -4L] # added jracine
    dfs <- c(0, aod$Df[-1L]) + object$rank; RSS <- aod$RSS
    n <- length(object$residuals)
    if(scale > 0) aic <- RSS/scale - n + k*dfs
    else aic <- n * log(RSS/n) + k*dfs
#    aod$CV <- aic
    o <- if(sorted) order(aod$CV) else seq_along(aod$CV)
    if(scale > 0) names(aod) <- c("Df", "Sum of Sq", "RSS", "Cp")
    test <- match.arg(test)
    if(test == "Chisq") {
        dev <- aod$"Sum of Sq"
        if(scale == 0) {
            dev <- n * log(RSS/n)
            dev <- dev[1L] - dev
            dev[1L] <- NA
        } else dev <- dev/scale
        df <- aod$Df
        nas <- !is.na(df)
        dev[nas] <- safe_pchisq(dev[nas], df[nas], lower.tail=FALSE)
        aod[, "Pr(Chi)"] <- dev
    } else if(test == "F") {
        rdf <- object$df.residual
        aod[, c("F Value", "Pr(F)")] <- Fstat(aod, aod$RSS[1L], rdf)
    }
    aod <- aod[o, ]
    head <- c("Single term additions", "\nModel:",
              deparse(as.vector(formula(object))))
    if(scale > 0)
        head <- c(head, paste("\nscale: ", format(scale), "\n"))
    class(aod) <- c("anova", "data.frame")
    attr(aod, "heading") <- head
    aod
}

addterm.negbin <- addterm.survreg <-
  function(object, ...)  addterm.default(object, ...)

addterm.glm <- function(object, scope, scale = 0, test = c("none", "Chisq", "F"),
                        k = 2, sorted = FALSE, trace = FALSE, ...)
{
  Fstat <- function(table, rdf) {
    dev <- table$Deviance
    df <- table$Df
    diff <- pmax(0, (dev[1L] - dev)/df)
    Fs <- (diff/df)/(dev/(rdf-df))
    Fs[df < .Machine$double.eps] <- NA
    P <- Fs
    nnas <- !is.na(Fs)
    P[nnas] <- safe_pf(Fs[nnas], df[nnas], rdf - df[nnas], lower.tail=FALSE)
    list(Fs=Fs, P=P)
  }
  if(missing(scope) || is.null(scope)) stop("no terms in scope")
  if(!is.character(scope))
    scope <- add.scope(object, update.formula(object, scope))
  if(!length(scope))
    stop("no terms in scope for adding to object")
  oTerms <- attr(terms(object), "term.labels")
  int <- attr(object$terms, "intercept")
  ns <- length(scope)
  dfs <- dev <- numeric(ns+1)
  names(dfs) <- names(dev) <- c("<none>", scope)
  add.rhs <- paste(scope, collapse = "+")
  add.rhs <- eval(parse(text = paste("~ . +", add.rhs)))
  new.form <- update.formula(object, add.rhs)
  oc <- object$call
  Terms <- terms(new.form)
  oc$formula <- Terms
  ## model.frame.glm looks at the terms part for the environment
  fob <- list(call = oc, terms=Terms)
  class(fob) <- class(object)
  x <- model.matrix(Terms, model.frame(fob, xlev = object$xlevels),
                    contrasts = object$contrasts)
  n <- nrow(x)
  oldn <- length(object$residuals)
  y <- object$y
  newn <- length(y)
  if(newn < oldn)
    warning(gettextf("using the %d/%d rows from a combined fit",
                     newn, oldn), domain = NA)
  wt <- object$prior.weights
  if(is.null(wt)) wt <- rep(1, n)
  Terms <- attr(Terms, "term.labels")
  asgn <- attr(x, "assign")
  ousex <- match(asgn, match(oTerms, Terms), 0L) > 0L
  if(int) ousex[1L] <- TRUE
  X <- x[, ousex, drop = FALSE]
  z <-  glm.fit(X, y, wt, offset=object$offset,
                family=object$family, control=object$control)
  dfs[1L] <- z$rank
  dev[1L] <- z$deviance
  ## workaround for PR#7842. terms.formula may have flipped interactions
  sTerms <- sapply(strsplit(Terms, ":", fixed=TRUE),
                   function(x) paste(sort(x), collapse=":"))
  for(tt in scope) {
    if(trace) {
      message("trying +", tt)
      utils::flush.console()
    }
    stt <- paste(sort(strsplit(tt, ":")[[1L]]), collapse=":")
    usex <- match(asgn, match(stt, sTerms), 0L) > 0L
    X <- x[, usex|ousex, drop = FALSE]
    z <-  glm.fit(X, y, wt, offset=object$offset,
                  family=object$family, control=object$control)
    dfs[tt] <- z$rank
    dev[tt] <- z$deviance
  }
  if (is.null(scale) || scale == 0)
    dispersion <- summary(object, dispersion = NULL)$dispersion
  else dispersion <- scale
  fam <- object$family$family
  if(fam == "gaussian") {
    if(scale > 0) loglik <- dev/scale - n
    else loglik <- n * log(dev/n)
  } else loglik <- dev/dispersion
  ## this cannot be correct jracine - but not reached ???
  aic <- loglik + k * dfs
  aic <- aic + (extractCV(object, k = k)[2L] - aic[1L]) # same baseline for CV
  dfs <- dfs - dfs[1L]
  dfs[1L] <- NA
  aod <- data.frame(Df = dfs, Deviance = dev, CV = aic,
                    row.names = names(dfs), check.names = FALSE)
  o <- if(sorted) order(aod$CV) else seq_along(aod$CV)
  if(all(is.na(aic))) aod <- aod[, -3]
  test <- match.arg(test)
  if(test == "Chisq") {
    dev <- pmax(0, loglik[1L] - loglik)
    dev[1L] <- NA
    LRT <- if(dispersion == 1) "LRT" else "scaled dev."
    aod[, LRT] <- dev
    nas <- !is.na(dev)
    dev[nas] <- safe_pchisq(dev[nas], aod$Df[nas], lower.tail=FALSE)
    aod[, "Pr(Chi)"] <- dev
  } else if(test == "F") {
    if(fam == "binomial" || fam == "poisson")
      warning(gettextf("F test assumes quasi%s family", fam),
              domain = NA)
    rdf <- object$df.residual
    aod[, c("F value", "Pr(F)")] <- Fstat(aod, rdf)
  }
  aod <- aod[o, ]
  head <- c("Single term additions", "\nModel:",
            deparse(as.vector(formula(object))))
  if(scale > 0)
    head <- c(head, paste("\nscale: ", format(scale), "\n"))
  class(aod) <- c("anova", "data.frame")
  attr(aod, "heading") <- head
  aod
}

addterm.mlm <- function(object, ...)
    stop("no addterm method implemented for \"mlm\" models")

dropterm <- function(object, ...) UseMethod("dropterm")

dropterm.default <-
  function(object, scope, scale = 0, test = c("none", "Chisq"),
           k = 2, sorted = FALSE, trace = FALSE, ...)
{
    tl <- attr(terms(object), "term.labels")
    if(missing(scope)) scope <- drop.scope(object)
    else {
        if(!is.character(scope))
            scope <- attr(terms(update.formula(object, scope)), "term.labels")
        if(!all(match(scope, tl, 0L)))
            stop("scope is not a subset of term labels")
    }
    ns <- length(scope)
    ans <- matrix(nrow = ns + 1L, ncol = 2L,
                  dimnames =  list(c("<none>", scope), c("df", "CV")))
    ans[1,  ] <- extractCV(object, scale, k = k, ...)
    env <- environment(formula(object))
    n0 <- length(object$residuals)
    for(i in seq(ns)) {
        tt <- scope[i]
        if(trace) {
      message("trying -", tt)
      utils::flush.console()
  }
        nfit <- update(object, as.formula(paste("~ . -", tt)),
                       evaluate = FALSE)
  nfit <- eval(nfit, envir=env) # was  eval.parent(nfit)
  ans[i+1, ] <- extractCV(nfit, scale, k = k, ...)
        if(length(nfit$residuals) != n0)
            stop("number of rows in use has changed: remove missing values?")
    }
    dfs <- ans[1L , 1L] - ans[, 1L]
    dfs[1L] <- NA
    aod <- data.frame(Df = dfs, CV = ans[,2])
    o <- if(sorted) order(aod$CV) else seq_along(aod$CV)
    test <- match.arg(test)
    if(test == "Chisq") {
        dev <- ans[, 2L] - k*ans[, 1L]
        dev <- dev - dev[1L] ; dev[1L] <- NA
        nas <- !is.na(dev)
        P <- dev
        P[nas] <- safe_pchisq(dev[nas], dfs[nas], lower.tail = FALSE)
        aod[, c("LRT", "Pr(Chi)")] <- list(dev, P)
    }
    aod <- aod[o, ]
    head <- c("Single term deletions", "\nModel:",
              deparse(as.vector(formula(object))))
    if(scale > 0)
        head <- c(head, paste("\nscale: ", format(scale), "\n"))
    class(aod) <- c("anova", "data.frame")
    attr(aod, "heading") <- head
    aod
}

dropterm.lm <-
  function(object, scope = drop.scope(object), scale = 0,
           test = c("none", "Chisq", "F"), k = 2, sorted = FALSE, ...)
{
#    aod <- stats:::drop1.lm(object, scope=scope, scale=scale)[, -4]
    aod <- drop1(object, scope=scope, scale=scale)[, -4]  
    dfs <-  object$rank - c(0, aod$Df[-1L]); RSS <- aod$RSS
    n <- length(object$residuals)
    aod$CV <- if(scale > 0)RSS/scale - n + k*dfs
    else n * log(RSS/n) + k*dfs
    o <- if(sorted) order(aod$CV) else seq_along(aod$CV)
    if(scale > 0) names(aod) <- c("Df", "Sum of Sq", "RSS", "Cp")
    test <- match.arg(test)
    if(test == "Chisq") {
        dev <- aod$"Sum of Sq"
        nas <- !is.na(dev)
        dev[nas] <- safe_pchisq(dev[nas]/scale, aod$Df[nas], lower.tail = FALSE)
        aod[, "Pr(Chi)"] <- dev
    } else if(test == "F") {
  dev <- aod$"Sum of Sq"
  dfs <- aod$Df
  rdf <- object$df.residual
  rms <- aod$RSS[1L]/rdf
  Fs <- (dev/dfs)/rms
  Fs[dfs < 1e-4] <- NA
  P <- Fs
  nas <- !is.na(Fs)
  P[nas] <- safe_pf(Fs[nas], dfs[nas], rdf, lower.tail=FALSE)
        aod[, c("F Value", "Pr(F)")] <- list(Fs, P)
    }
    aod <- aod[o, ]
    head <- c("Single term deletions", "\nModel:",
              deparse(as.vector(formula(object))))
    if(scale > 0)
        head <- c(head, paste("\nscale: ", format(scale), "\n"))
    class(aod) <- c("anova", "data.frame")
    attr(aod, "heading") <- head
    aod
}

dropterm.mlm <- function(object, ...)
  stop("dropterm not implemented for \"mlm\" fits")

dropterm.glm <-
  function(object, scope, scale = 0, test = c("none", "Chisq", "F"),
           k = 2, sorted = FALSE, trace = FALSE, ...)
{
    x <- model.matrix(object)
    n <- nrow(x)
    asgn <- attr(x, "assign")
    tl <- attr(object$terms, "term.labels")
    if(missing(scope)) scope <- drop.scope(object)
    else {
        if(!is.character(scope))
            scope <- attr(terms(update.formula(object, scope)), "term.labels")
        if(!all(match(scope, tl, 0L)))
            stop("scope is not a subset of term labels")
  }
    ns <- length(scope)
    ndrop <- match(scope, tl)
    rdf <- object$df.residual
    chisq <- object$deviance
    dfs <- numeric(ns)
    dev <- numeric(ns)
    y <- object$y
    if(is.null(y)) {
        y <- model.response(model.frame(object))
        if(!is.factor(y)) storage.mode(y) <- "double"
    }
    wt <- object$prior.weights
    if(is.null(wt)) wt <- rep.int(1, n)
    for(i in 1L:ns) {
        if(trace) {
      message("trying -", scope[i])
      utils::flush.console()
  }
        ii <- seq_along(asgn)[asgn == ndrop[i]]
        jj <- setdiff(seq(ncol(x)), ii)
        z <-  glm.fit(x[, jj, drop = FALSE], y, wt, offset=object$offset,
                      family=object$family, control=object$control)
        dfs[i] <- z$rank
        dev[i] <- z$deviance
    }
    scope <- c("<none>", scope)
    dfs <- c(object$rank, dfs)
    dev <- c(chisq, dev)
    dispersion <- if (is.null(scale) || scale == 0)
  summary(object, dispersion = NULL)$dispersion
    else scale
    fam <- object$family$family
    loglik <-
        if(fam == "gaussian") {
            if(scale > 0) dev/scale - n else n * log(dev/n)
        } else dev/dispersion
    aic <- loglik + k * dfs
    dfs <- dfs[1L] - dfs
    dfs[1L] <- NA
    aic <- aic + (extractCV(object, k = k)[2L] - aic[1L])
    aod <- data.frame(Df = dfs, Deviance = dev, CV = aic,
                      row.names = scope, check.names = FALSE)
    o <- if(sorted) order(aod$CV) else seq_along(aod$CV)
    if(all(is.na(aic))) aod <- aod[, -3]
    test <- match.arg(test)
    if(test == "Chisq") {
        dev <- pmax(0, loglik - loglik[1L])
        dev[1L] <- NA
        nas <- !is.na(dev)
        LRT <- if(dispersion == 1) "LRT" else "scaled dev."
        aod[, LRT] <- dev
        dev[nas] <- safe_pchisq(dev[nas], aod$Df[nas], lower.tail=FALSE)
        aod[, "Pr(Chi)"] <- dev
    } else if(test == "F") {
        if(fam == "binomial" || fam == "poisson")
            warning(gettextf("F test assumes 'quasi%s' family", fam),
                    domain = NA)
  dev <- aod$Deviance
  rms <- dev[1L]/rdf
        dev <- pmax(0, dev - dev[1L])
  dfs <- aod$Df
  rdf <- object$df.residual
  Fs <- (dev/dfs)/rms
  Fs[dfs < 1e-4] <- NA
  P <- Fs
  nas <- !is.na(Fs)
  P[nas] <- safe_pf(Fs[nas], dfs[nas], rdf, lower.tail=FALSE)
  aod[, c("F value", "Pr(F)")] <- list(Fs, P)
    }
    aod <- aod[o, ]
    head <- c("Single term deletions", "\nModel:",
              deparse(as.vector(formula(object))))
    if(scale > 0)
        head <- c(head, paste("\nscale: ", format(scale), "\n"))
    class(aod) <- c("anova", "data.frame")
    attr(aod, "heading") <- head
    aod
}

dropterm.negbin <- dropterm.survreg <-
    function(object, ...) dropterm.default(object, ...)
