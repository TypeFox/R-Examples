cox.kmi <- function(formula, imp.data, df.complete = Inf, ...) {
    if (!inherits(imp.data, "kmi")) {
        stop("'imp.data' must be of class 'kmi'")
    }
    call <- match.call()
    info <- imp.data$info # that's where we have the column names (time, event)
    result <- lapply(seq_along(imp.data$imputed.data), function(i) {
        daten <- imp.data$original.data
        daten[, info[1]] <- imp.data$imputed.data[[i]][, 1]
        daten[, info[2]] <- imp.data$imputed.data[[i]][, 2]
        tmp <- coxph(formula, data = daten, ...)
        tmp
    })
    res <- MIcombine(result, df.complete = df.complete) ## that's a nice function
    zzz <- list(coefficients = res$coefficients,
                variance = res$variance,
                nimp = res$nimp,
                df = res$df,
                call = call,
                individual.fit = result)
    class(zzz) <- "cox.kmi"
    zzz
}
