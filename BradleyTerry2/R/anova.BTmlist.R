anova.BTmlist <- function (object, ..., dispersion = NULL, test = NULL) {
    ## Pass on if no random effects
    fixed <- unlist(lapply(object, function(x) is.null(x$random)))
    if (!all(!fixed))
        stop("Models must have the same random effects structure")

    responses <- as.character(lapply(object, function(x) {
        deparse(formula(terms(x))[[2]])
    }))
    sameresp <- responses == responses[1]
    if (!all(sameresp)) {
        object <- object[sameresp]
        warning("models with response ", deparse(responses[!sameresp]),
                " removed because response differs from model 1")
    }
    ns <- sapply(object, function(x) length(fitted(x)))
    if (any(ns != ns[1]))
        stop("models were not all fitted to the same size of dataset")
    nmodels <- length(object)

    ncoefs <- sapply(object, function(x) length(na.omit(coef(x)))) #omit aliased
    labels <- lapply(object, function(x) x$term.labels)
    stat <- numeric(nmodels)
    for (i in 2:nmodels) {
        descending <- ncoefs[i] < ncoefs[i - 1]
        bigger <- i - descending
        smaller <- i - !descending
        if (!all(labels[[smaller]] %in% labels[[bigger]]))
            stop("models are not nested")
        ind <- !(labels[[bigger]] %in% labels[[smaller]])
        stat[i] <- t(coef(object[[bigger]])[ind]) %*%
            chol2inv(chol(vcov(object[[bigger]], dispersion = dispersion)[ind, ind])) %*%
                coef(object[[bigger]])[ind] #vcov should deal with dispersion != 1
    }
    stat[1] <- NA
    table <- data.frame(stat, c(NA, diff(ncoefs)))
    variables <- lapply(object, function(x) paste(deparse(formula(x)),
                                                  collapse = "\n"))
    dimnames(table) <- list(1:nmodels, c("Statistic", "Df"))
    title <- paste("Sequential Wald Tests\n\n",
                   "Response: ", responses[1], "\n", sep = "")
    topnote <- paste("Model ", format(1:nmodels), ": ", variables,
                     sep = "", collapse = "\n")
    if (!is.null(test)) {
        ## Assume dispersion fixed at one - if dispersion estimated, would use
        ## "residual" df from larger model in each comparison
        df.dispersion <- Inf
        if (test == "F" && df.dispersion == Inf) {
            fam <- object[[1]]$family$family
            if (fam == "binomial" || fam == "poisson")
                warning(gettextf("using F test with a '%s' family is inappropriate",
                  fam), domain = NA, call. = FALSE)
            else warning("using F test with a fixed dispersion is inappropriate")
        }
        table <- switch(test, Chisq = {
            dfs <- table[, "Df"]
            vals <- table[, "Statistic"]
            vals[dfs %in% 0] <- NA
            cbind(table, `P(>|Chi|)` = pchisq(vals, abs(dfs), lower.tail = FALSE))
        }, F = {
            dfs <- table[, "Df"]
            Fvalue <- table[, "Statistic"]/abs(dfs)
            Fvalue[dfs %in% 0] <- NA
            cbind(table, F = Fvalue, `Pr(>F)` =
                  pf(Fvalue, abs(dfs), df.dispersion, lower.tail = FALSE))
        })
    }
    structure(table, heading = c(title, topnote), class = c("anova",
        "data.frame"))
}





