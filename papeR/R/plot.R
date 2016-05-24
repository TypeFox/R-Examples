## plot function for (labeled) data frames

## now specify modified plot function for data frames
plot.ldf <- function(x, variables = names(x),
                     labels = TRUE, by = NULL,
                     with = NULL, regression.line = TRUE,
                     line.col = "red", ...) {

    if (is.numeric(variables)) {
        variables <- names(x)[variables]
    }

    if (!is.null(with)) {
        if (!is.null(by))
            stop("One can only specify either ", sQuote("by"), " or ",
                 sQuote("with"))
        by <- with
    }

    if (is.numeric(by)) {
        by <- names(x)[by]
    }

    if (!all(c(by, variables) %in% names(x)))
        stop("(Some of the) specified variables are not available in x")

    ## set up labels
    if (is.null(labels)) {
        labels <- variables
    } else {
        if (is.logical(labels) && labels) {
            labels <- labels(x, which = variables)
            if (!is.null(by))
                grp_label <- labels(x, which = by)
        } else {
            if (length(variables) != length(labels))
                stop(sQuote("variables"), " and ", sQuote("labels"),
                     " must have the same length")
        }
    }

    if (!is.null(by)) {
        if(!is.factor(x[, by]) && !is.numeric(x[, by]))
            stop(sQuote("by"), " must specify a factor or numeric variable")
        if (by %in% variables) {
            idx <- variables != by
            variables <- variables[idx]
            labels <- labels[idx]
        }
        by_var <- x[, by]
    }

    x <- x[, variables, drop = FALSE]

    ## get numerical variables
    num <- mySapply(x, is.numeric)
    fac <- mySapply(x, is.factor)

    ## if anything else is present (not num or fac)
    if (!all(num | fac))
        warning("Only numeric or factor variables are plotted")

    which.num <- which(num)
    which.fac <- which(fac)

    if (is.null(by)) {
        for (i in which.num) {
            boxplot(x[, i], main = variables[i], ylab = labels[i], ...)
        }
        for (i in which.fac) {
            barplot(table(x[, i]),
                    main = variables[i], ylab = labels[i], ...)
        }
    } else {
        grp_label <- ifelse(!is.null(grp_label), grp_label, by)
        if (is.factor(by_var)) {
            for (i in which.num) {
                cc <- complete.cases(x[, i], by_var)
                tmp_by_var <- by_var[cc, drop = TRUE]
                boxplot(x[cc, i] ~ tmp_by_var, main = variables[i],
                        ylab = labels[i], xlab = grp_label, ...)
            }
            for (i in which.fac) {
                cc <- complete.cases(x[, i], by_var)
                tmp_by_var <- by_var[cc, drop = TRUE]
                plot(tmp_by_var, x[cc, i], main = variables[i],
                     ylab = labels[i], xlab = grp_label, ...)
            }
        } else {  ## i.e. is.numeric(by_var)
            for (i in which.num) {
                cc <- complete.cases(x[, i], by_var)
                tmp_by_var <- by_var[cc, drop = TRUE]
                graphics::plot.default(x[cc, i], tmp_by_var, main = variables[i],
                                       xlab = labels[i], ylab = grp_label, ...)
                if (regression.line)
                    abline(lm(tmp_by_var ~ x[cc, i]), col = line.col)
            }
            for (i in which.fac) {
                cc <- complete.cases(x[, i], by_var)
                tmp_by_var <- by_var[cc, drop = TRUE]
                boxplot(tmp_by_var ~ x[cc, i], main = variables[i],
                        xlab = labels[i], ylab = grp_label, ...)
            }
        }
    }
}
