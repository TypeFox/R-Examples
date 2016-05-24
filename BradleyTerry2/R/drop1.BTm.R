drop1.BTm <- function(object, scope, scale = 0, test = c("none", "Chisq", "F"),
                      ...) {
    x <- model.matrix(object)

    ## Pass on if no random effects
    if (is.null(object$random)){
        object$x <- x
        attr(object$x, "assign") <- object$assign
        object$terms <- terms(object$formula)
        return(NextMethod())
    }

    form <- formula(object)

    if (missing(scope))
        scope <- drop.scope(nobars(form))
     else {
        if (!is.character(scope)) {
            srandom <- findbars(scope[[2]])
            if (length(srandom))
                stop("Scope should not include random effects.")
            scope <- attr(terms(update.formula(form, scope)),
                          "term.labels")
        }
        if (!all(match(scope, terms(form), 0L) > 0L))
            stop("scope is not a subset of term labels")
    }

    asgn <- object$assign

    coefs <- coef(object)
    if (scale == 0) dispersion <- 1
    else dispersion <- scale
    vc <- vcov(object, dispersion = dispersion) #vcov should deal with dispersion != 1

    sTerms <- sapply(strsplit(scope, ":", fixed = TRUE),
                     function(x) paste(sort(x), collapse = ":"))
    stat <- df <- numeric(length(scope))
    names(stat) <- names(df) <- as.character(sapply(scope, as.name))
    tryerror <- FALSE
    for (i in seq(scope)) {
        stt <- paste(sort(strsplit(scope[i], ":")[[1]]), collapse = ":")
        usex <- match(asgn, match(stt, sTerms), 0) > 0
        trystat <- try(t(coefs[usex]) %*% chol2inv(chol(vc[usex, usex])) %*%
                       coefs[usex], silent = TRUE)
        if (inherits(trystat, "try-error")) {
            stat[i] <- df[i] <- NA
            tryerror <- TRUE
        }
        else {
            stat[i] <- trystat
            df[i] <- sum(usex)
        }
    }
    table <- data.frame(stat, df)
    dimnames(table) <- list(names(df), c("Statistic", "Df"))
    title <- "Single term deletions\n"
    topnote <- gsub("\\s+", " ", paste("Model: ",
                                     paste(deparse(as.vector(formula(object))),
                                           collapse = ""),
                     if (scale > 0) paste("\nscale: ", format(scale), "\n"),
                     if (tryerror)
                     "\n\nTest statistic unestimable for at least one term"),
                    perl = TRUE)
    test <- match.arg(test)
    if (test == "Chisq") {
        dfs <- table[, "Df"]
        vals <- table[, "Statistic"]
        vals[dfs %in% 0] <- NA
        table <- cbind(table, `P(>|Chi|)` = pchisq(vals, abs(dfs),
                              lower.tail = FALSE))
    }
    else if (test == "F") {
        ## Assume dispersion fixed at one - if dispersion estimated, would use
        ## "residual" df from larger model in each comparison
        df.dispersion <- Inf
        if (df.dispersion == Inf) {
            fam <- object[[1]]$family$family
            if (fam == "binomial" || fam == "poisson")
                warning(gettextf("using F test with a '%s' family is inappropriate",
                                 fam), domain = NA, call. = FALSE)
            else warning("using F test with a fixed dispersion is inappropriate")
        }
        dfs <- table[, "Df"]
        Fvalue <- table[, "Statistic"]/abs(dfs)
        Fvalue[dfs %in% 0] <- NA
        table <- cbind(table, F = Fvalue, `Pr(>F)` =
                       pf(Fvalue, abs(dfs), df.dispersion,
                          lower.tail = FALSE))
    }
    structure(table, heading = c(title, topnote), class = c("anova",
        "data.frame"))
}
