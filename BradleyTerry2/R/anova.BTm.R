anova.BTm <- function (object, ..., dispersion = NULL, test = NULL)
{
    ## Only list models in ...
    dotargs <- list(...)
    named <- if (is.null(names(dotargs)))
        rep(FALSE, length(dotargs))
    else (names(dotargs) != "")
    if (any(named))
        warning("the following arguments to 'anova.BTm' are invalid and dropped: ",
            paste(deparse(dotargs[named]), collapse = ", "))
    dotargs <- dotargs[!named]
    is.BTm <- unlist(lapply(dotargs, function(x) inherits(x, "BTm")))
    dotargs <- dotargs[is.BTm]

    ## Compare list of models
    models <- c(list(object), dotargs)
    if (length(dotargs) > 0){
        fixed <- unlist(lapply(models, function(x) is.null(x$random)))
        if (all(fixed)) {
            variables <- lapply(models, function(x) paste(deparse(formula(x)),
                                                          collapse = "\n"))
            models <- lapply(models, function(x) {
                x$formula <- formula(x$terms)
                class(x) <- setdiff(class(x), "BTm")
                x})
            call <- match.call()
            anova.table <- do.call("anova", c(models, call$dispersion, call$test))
            attr(anova.table, "heading") <-
                c(paste("Analysis of Deviance Table\n\n",
                        "Response: ", deparse(object$call$outcome, 500), "\n", sep = ""),
                  paste("Model ", format(seq(models)), ": ", variables,
                        sep = "", collapse = "\n"))
            return(anova.table)
        }
        else
            return(anova.BTmlist(c(list(object), dotargs),
                                 dispersion = dispersion,  test = test))
    }

    X <- model.matrix(object)
    Z <- object$random

    sep <- 0 %in% object$assign

    ## Passing on to glm when no random effects
    if (is.null(Z)) {
        object$x <- X
        attr(object$x, "assign") <- object$assign + sep
        attr(object$terms, "term.labels") <- c("[sep]"[sep], object$term.labels)
        anova.table <- NextMethod()
        attr(anova.table, "heading") <-
            paste("Analysis of Deviance Table", "\n\nModel: ",
                  object$family$family, ", link: ", object$family$link,
                  "\n\nResponse: ", deparse(object$call$outcome, 500),
                  "\n\nTerms added sequentially (first to last)\n\n",
                  sep = "")
        if (sep) {
            anova.table <- anova.table[-1,]
            rownames(anova.table)[1] <- "NULL"
            anova.table[1, 1:2] <- NA
        }
        return(anova.table)
    }

    varseq <- object$assign

    nvars <- max(0, varseq)
    stat <- df <- numeric(nvars)
    tryerror <- FALSE
    if (nvars > 1) {
        y <- object$y
        ## Extension to further methods
        method <- object$method
        if (!is.function(method))
            method <- get(method, mode = "function")
        control <- object$control
        control$trace <- FALSE
        for (i in 1:(nvars - 1)) {
            fit <- method(X = X[, varseq <= i, drop = FALSE], y = y, Z = Z,
                          weights = object$prior.weights, start = object$start,
                          offset = object$offset, family = object$family,
                          control = control,
                          sigma = object$call$sigma,
                          sigma.fixed = object$sigma.fixed)
            class(fit) <- oldClass(object)
            ind <- (varseq == i)[varseq <= i]
            trystat <- try(t(coef(fit)[ind]) %*%
                           chol2inv(chol(suppressMessages(vcov(fit, dispersion = dispersion))[ind, ind])) %*%
                           coef(fit)[ind], silent = TRUE) #vcov should deal with dispersion != 1
            if (inherits(trystat, "try-error")) {
                stat[i] <- df[i] <- NA
                tryerror <- TRUE
            }
            else {
                stat[i] <- trystat
                df[i] <- sum(ind)
            }
        }
    }
    ind <- varseq == nvars
    trystat <- try(t(coef(object)[ind]) %*% chol2inv(chol(object$varFix[ind, ind])) %*%
                   coef(object)[ind], silent = TRUE)
    if (inherits(trystat, "try-error")) {
        stat[nvars] <- df[nvars] <- NA
        tryerror <- TRUE
    }
    else {
        stat[nvars] <- trystat
        df[nvars] <- sum(ind)
    }
    table <- data.frame(c(NA, stat), c(NA, df))
    dimnames(table) <- list(c("NULL", object$term.labels), c("Statistic", "Df"))
    title <- paste("Sequential Wald Tests", "\n\nModel: ",
                   object$family$family, ", link: ", object$family$link,
                   "\n\nResponse: ", deparse(object$call$outcome, 500),
                   "\n\nPredictor: ", paste(formula(object), collapse = ""),
                   "\n\nTerms added sequentially (first to last)",
                   if (tryerror)
                   "\n\nTest statistic unestimable for at least one term",
                   "\n", sep = "")

    ## Assume dispersion fixed at one - if dispersion estimated, would use
    ## "residual" df from larger model in each comparison
    df.dispersion <- Inf
    if (!is.null(test)) {
        if (test == "F" && df.dispersion == Inf) {
            fam <- object$family$family
            if (fam == "binomial" || fam == "poisson")
                warning(gettextf("using F test with a %s family is inappropriate",
                  fam), domain = NA)
            else warning("using F test with a fixed dispersion is inappropriate")
        }
        table <- switch(test, Chisq = {
            dfs <- table[, "Df"]
            vals <- table[, "Statistic"]
            vals[dfs %in% 0] <- NA
            cbind(table, `P(>|Chi|)` = pchisq(vals, dfs, lower.tail = FALSE))
        }, F = {
            dfs <- table[, "Df"]
            Fvalue <- table[, "Statistic"]/dfs
            Fvalue[dfs %in% 0] <- NA
            cbind(table, F = Fvalue, `Pr(>F)` =
                  pf(Fvalue, dfs, df.dispersion, lower.tail = FALSE))
        })
    }
    structure(table, heading = title, class = c("anova", "data.frame"))
}
