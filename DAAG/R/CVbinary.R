CVbinary <-
function (obj, rand = NULL, nfolds = 10, print.details = TRUE)
{
    data <- obj$data
    m <- dim(data)[1]
    if (is.null(rand))
        rand <- sample(nfolds, m, replace = TRUE)
    form <- formula(obj)
    yvar <- all.vars(form)[1]
    obs <- data[, yvar]
    ival <- unique(rand)
    fam <- obj$family$family
    hat <- predict(glm(form, data, family = fam), type = "response")
    cvhat <- rep(0, length(rand))
    if (print.details)
        cat("\nFold: ")
    for (i in ival) {
        if (print.details)
            cat("", i)
        if (i%%20 == 0)
            cat("\n")
        here <- i != rand
        i.glm <- glm(form, data = data[here, ], family = fam)
        cvhat[!here] <- predict(i.glm, newdata = data[!here,
            ], family = fam, type = "response")
    }
    if (is.factor(obs)) {
        lev <- levels(obs)
        hat <- lev[round(hat) + 1]
        cvhat <- lev[round(cvhat) + 1]
        acc.internal <- sum(obs == hat)/m
        acc.cv <- sum(obs == cvhat)/m
    }
    else {
        acc.internal <- sum(obs == round(hat))/m
        acc.cv <- sum(obs == round(cvhat))/m
    }
    if (print.details) {
        cat("\nInternal estimate of accuracy =", round(acc.internal,
            3))
        cat("\nCross-validation estimate of accuracy =", round(acc.cv,
            3))
        cat("\n")
    }
    invisible(list(cvhat = cvhat, internal = hat, training=hat, acc.cv = acc.cv,
                   acc.internal = acc.internal, acc.training=acc.internal))
}
