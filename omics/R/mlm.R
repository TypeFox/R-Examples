mlm <- function(formula, data, vars, save.residuals=FALSE) {
    if (missing(data)) {
        data <- NULL
    }
    formula <- as.formula(terms(formula, data=data))

    Y <- get(response.name(formula), envir=environment(formula))

    formula <- update(formula, "NULL ~ .")
    mf <- model.frame(formula, data, na.action=na.pass, drop.unused.levels=TRUE)
    mm <- model.matrix(formula, mf)

    labs <- labels(terms(formula))
    if (missing(vars)) {
        vars <- labs
    }

    idx <- which(attr(mm, "assign") %in% match(vars, labs))
    if (length(idx) == 0 & !save.residuals) {
        stop("No variables selected")
    }
    save.coefs <- !(length(idx) == 0 & save.residuals)
    vars <- colnames(mm)[idx]
    colnames(mm) <- sprintf("V%d", 1:ncol(mm))
    new.vars <- colnames(mm)[idx]
    mm <- as.data.frame(mm)
    formula <- as.formula(sprintf("y ~ %s - 1",
        paste0(colnames(mm), collapse=" + ")
    ))

    ns <- rep(NA, ncol(Y))

    if (save.coefs) {
        coefs <- array(NA, c(ncol(Y), length(vars), 3),
                       dimnames=list(colnames(Y), vars,
                                     c("coef", "coef.se", "pval")))
    }

    if (save.residuals) {
        residuals <- matrix(NA, nrow(Y), ncol(Y), dimnames=dimnames(Y))
    }

    for (i in 1:ncol(Y)) {
        mm$y <- Y[,i]

        model <- try(lm(formula, data=mm, na.action=na.exclude), silent=TRUE)
        if (inherits(model, "try-error")) {
            next
        }

        ns[i] <- nobs(model)

        if (save.coefs) {
            tmp <- try(coef(summary(model)), silent=TRUE)
            if (inherits(tmp, "try-error")) {
                next
            }

            for (j in 1:length(vars)) if (new.vars[j] %in% rownames(tmp)) {
                coefs[i,vars[j],"coef"] <- tmp[new.vars[j],"Estimate"]
                coefs[i,vars[j],"coef.se"] <- tmp[new.vars[j],"Std. Error"]
                coefs[i,vars[j],"pval"] <- tmp[new.vars[j],"Pr(>|t|)"]
            }
        }

        if (save.residuals) {
            residuals[,i] <- resid(model, type="response")
        }
    }

    if (save.coefs & length(vars) == 1) {
        coefs <- as.data.frame(coefs[,1,])
    }

    tmp <- list(
        nobs=ns
    )
    if (save.coefs) {
        tmp$coefficients <- coefs
    }
    if (save.residuals) {
        tmp$residuals <- residuals
    }
    tmp
}
