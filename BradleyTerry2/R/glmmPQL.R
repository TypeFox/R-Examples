glmmPQL <- function(fixed, random = NULL, family = binomial, data = NULL,
                    subset = NULL, weights = NULL, offset = NULL,
                    na.action = NULL,  start = NULL, etastart = NULL,
                    mustart = NULL, control = glmmPQL.control(...),
                    sigma = 0.1, sigma.fixed = FALSE, model = TRUE,
                    x = FALSE, contrasts = NULL, ...) {
    call <-  match.call()
    nm <- names(call)[-1]

    if (is.null(random)) {
        keep <- is.element(nm, c("family", "data", "subset", "weights",
                                 "offset", "na.action"))
        for (i in nm[!keep]) call[[i]] <- NULL
        call$formula <- fixed
        environment(call$formula) <- environment(fixed)
        call[[1]] <- as.name("glm")
        return(eval.parent(call))
    }

    modelTerms <- terms(fixed, data = data)
    modelCall <- as.list(match.call(expand.dots = FALSE))
    argPos <- match(c("data", "subset", "na.action", "weights", "offset"),
                    names(modelCall), 0)
    modelData <- as.call(c(model.frame, list(formula = modelTerms,
                                             drop.unused.levels = TRUE),
                           modelCall[argPos]))
    modelData <- eval(modelData, parent.frame())

    if (!is.null(modelCall$subset))
        Z <- random[eval(modelCall$subset, data, parent.frame()),]
    else Z <- random

    if (!is.null(attr(modelData, "na.action")))
        Z <- Z[-attr(modelData, "na.action"),]

    nObs <- nrow(modelData)
    y <- model.response(modelData, "numeric")
    if (is.null(y))
        y <- rep(0, nObs)
    weights <- as.vector(model.weights(modelData))
    if (!is.null(weights) && any(weights < 0))
        stop("negative weights are not allowed")
    if (is.null(weights))
        weights <- rep.int(1, nObs)
    offset <- as.vector(model.offset(modelData))
    if (is.null(offset))
        offset <- rep.int(0, nObs)

    if (is.character(family))
        family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family))
        family <- family()
    if (is.null(family$family)) {
        print(family)
        stop("`family' not recognized")
    }
    if (family$family == "binomial") {
        if (is.factor(y) && NCOL(y) == 1)
            y <- y != levels(y)[1]
        else if (NCOL(y) == 2) {
            n <- y[, 1] + y[, 2]
            y <- ifelse(n == 0, 0, y[, 1]/n)
            weights <- weights * n
        }
    }

    ## Use GLM to estimate fixed effects
    empty <- is.empty.model(modelTerms)
    if (!empty)
        X <- model.matrix(fixed, data = modelData, contrasts)
    else
        X <- matrix(, nObs, 0)
    fit <- glmmPQL.fit(X = X, y = y, Z = Z, weights = weights, start = start,
                       etastart = etastart, mustart = mustart, offset = offset,
                       family = family, control = control, sigma = sigma,
                       sigma.fixed = sigma.fixed, ...)
    if (sum(offset) && attr(modelTerms, "intercept") > 0) {
        fit$null.deviance <- glm.fit(x = X[, "(Intercept)", drop = FALSE],
            y = y, weights = weights, offset = offset, family = family,
            control = control, intercept = TRUE)$deviance
    }
    if (model)
        fit$model <- modelData
    fit$na.action <- attr(modelData, "na.action")
    if (x)
        fit$x <- X
    fit <- c(fit, list(call = call, formula = fixed, random = random,
                       terms = modelTerms,
                       data = data, offset = offset, control = control,
                       method = "glmmPQL.fit", contrasts = attr(X, "contrasts"),
                       xlevels = .getXlevels(modelTerms, modelData)))
    class(fit) <- c("BTglmmPQL", "glm", "lm")
    fit
}

