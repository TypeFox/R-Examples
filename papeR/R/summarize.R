################################################################################
##  Author: Benjamin Hofner, benjamin.hofner@fau.de

## define summarize
summarize <- summarise <- function(data, type = c("numeric", "factor"),
    variables = names(data), variable.labels = labels, labels = NULL, group = NULL,
    test = !is.null(group), colnames = NULL, digits = NULL, digits.pval = 3,
    smallest.pval = 0.001, sep = NULL, sanitize = TRUE, drop = TRUE,
    show.NAs = any(is.na(data[, variables])), ...) {

    type <- match.arg(type)

    ## get call
    cll <- match.call()
    ## modify call
    cll[[1]] <- as.name(paste0("summarize_", type))
    cll$type <- NULL

    if (is.null(digits)) {
        cll$digits <- NULL
    }
    if (is.null(sep)) {
        cll$sep <- NULL
    }
    ## and evaluate the modified call
    eval(cll)
}

################################################################################
### for backward compatibility define
latex.table.cont <- function(...,
                             caption = NULL, label = NULL,
                             table = c("tabular", "longtable"), align = NULL,
                             floating = FALSE, center = TRUE) {

    tab <- summarize_numeric(...)
    table <- match.arg(table)
    print(xtable(tab, caption = caption, label = label, align = align),
          floating = floating, latex.environments = ifelse(center, "center", c()),
          tabular.environment = table)
    message("  This function exists for backward compatibility.\n  Consider using ",
            sQuote('xtable(summarize(..., type = "numeric"))'), " instead.")
}

latex.table.fac <- function(...,
                            caption = NULL, label = NULL,
                            table = c("tabular", "longtable"), align = NULL,
                            floating = FALSE, center = TRUE) {

    tab <- summarize_factor(...)
    table <- match.arg(table)
    print(xtable(tab, caption = caption, label = label, align = align),
          floating = floating, latex.environments = ifelse(center, "center", c()),
          tabular.environment = table)
    message("  This function exists for backward compatibility.\n  Consider using ",
            sQuote('xtable(summarize(..., type = "factor"))'), " instead.")

}
################################################################################


################################################################################
# LaTeX Tables with Descriptves for Continuous Variables
summarize_numeric <- function(data, variables = names(data),
                              variable.labels = labels, labels = NULL, group = NULL,
                              test = !is.null(group), colnames = NULL,
                              digits = 2, digits.pval = 3,
                              smallest.pval = 0.001, sep = !is.null(group),
                              sanitize = TRUE, drop = TRUE,
                              show.NAs = any(is.na(data[, variables])),
                              count = TRUE, mean_sd = TRUE,
                              quantiles = TRUE, incl_outliers = TRUE, ...) {

    if (is.null(variable.labels)) {
        variable.labels <- variables
    } else {
        if (is.logical(variable.labels) && variable.labels) {
            variable.labels <- labels(data, which = variables)
        } else {
            if (length(variables) != length(variable.labels))
                stop(sQuote("variables"), " and ", sQuote("variable.labels"),
                     " must have the same length")
        }
    }

    if (!is.null(group)) {
        if(!is.factor(data[, group]))
            stop(sQuote("group"), " must be a factor variable")
        if (group %in% variables) {
            idx <- variables != group
            variables <- variables[idx]
            variable.labels <- variable.labels[idx]
        }
        ## remove observations with missing group:
        if (any(is.na(data[, group]))) {
            warning("Removed observations with missing group")
            data <- data[!is.na(data[, group]), ]
        }
        group_var <- data[, group]
    }
    ## get numerical variables
    num <- mySapply(data[, variables], function(x)
        is.numeric(x) || inherits(x, "Date"))
    date <- mySapply(data[, variables], function(x)
        inherits(x, "Date"))

    ## drop missings
    if (drop) {
        compl.missing <- mySapply(data[, variables], function(x) all(is.na(x)))
        num <- num & !compl.missing
        date <- date & !compl.missing
    }

    ## if not any is TRUE (i.e. all are FALSE):
    if (!any(count, mean_sd, quantiles))
        stop("Nothing to compute. All quantities are set to FALSE.")
    ## if all variables are factors:
    if (all(!num))
        stop("None of the variables is numeric or all variables are missing. Nothing to compute.")
    ## if factors are dropped:
    if (any(!num))
        message("Factors are dropped from the summary")

    ## subset variables to non-factors only
    variables <- variables[num]
    variable.labels <- variable.labels[num]

    ## setup results object
    sums <- data.frame(variable = variable.labels, group = NA, blank = "",
                       N=NA, Missing = NA, blank_1 = "",
                       Mean=NA, SD=NA, blank_2 = "",
                       Min=NA, Q1=NA, Median=NA, Q3=NA, Max=NA, var = variables,
                       stringsAsFactors = FALSE)

    if (!is.null(group)) {
        sums <- sums[rep(1:nrow(sums), each = nlevels(group_var)), ]
        sums$group <- levels(group_var)
    } else {
        ## drop group variable
        sums$group <- NULL
        sums$blank <- NULL
    }

    ## compute statistics
    for (i in 1:nrow(sums)) {
        if (!is.null(group)) {
            idx <- c(4:5, 7:8, 10:14)
            sums[i, idx] <- compute_summary(data[, sums$var[i]],
                                            group_var = group_var,
                                            group = sums$group[i],
                                            incl_outliers = incl_outliers,
                                            digits = digits)
        } else {
            idx <- c(2:3, 5:6, 8:12)
            sums[i, idx] <- compute_summary(data[, sums$var[i]],
                                            incl_outliers = incl_outliers,
                                            digits = digits)
        }
    }

    if (!is.null(group)) {
        if (!is.character(test) && test)
            test <- "t.test"

        if (all(is.character(test))) {
            if (length(test) == 1)
                test <- rep(test, length(variables))

            pval <- rep(NA, length(variables))
            for (i in 1:length(variables)) {
                fm <- as.formula(paste(variables[i], " ~ ", group))
                pval[i] <- do.call(test[i], list(formula = fm, data = data))$p.value
            }
            ## make sure rounding is to digits.pval digits
            pval <- format.pval(round(pval, digits = digits.pval),
                                eps = smallest.pval)
            ## make sure not to drop trailing zeros
            pval2 <- suppressWarnings(as.numeric(pval))
            pval[!is.na(pval2)] <- sprintf(paste0("%0.", digits.pval, "f"),
                                           pval2[!is.na(pval2)])
            pval[is.na(pval)] <- ""
            sums$blank_p <- ""
            sums$p.value <- ""
            sums$p.value[!duplicated(sums$var)] <- pval
        }
    }

    ## remove superfluous variables
    sums$var <- NULL
    if (!show.NAs) {
        sums$Missing <- NULL
    }
    if (!is.null(group)) {
        names(sums)[names(sums) == "group"] <- labels(data, group)
    }

    sums <- set_options(sums, sep = sep, sanitize = sanitize,
                        count = count, mean_sd = mean_sd, quantiles = quantiles,
                        colnames = colnames, class = "summarize.numeric")
    prettify(sums)
}

compute_summary <- function(data, ...)
    UseMethod("compute_summary")

compute_summary.default <- function(data, group_var = NULL, group = NULL,
                                    incl_outliers, digits) {
    if (!is.null(group)) {
        data <- data[group_var == group]
    }

    res <- data.frame(N=NA, Missing = NA, Mean=NA, SD=NA,
                      Min=NA, Q1=NA, Median=NA, Q3=NA, Max=NA)

    res["N"] <- sum(!is.na(data))
    res["Missing"] <- sum(is.na(data))
    res["Mean"] <- round(mean(data, na.rm=TRUE), digits = digits)
    res["SD"] <- round(sd(data, na.rm=TRUE), digits = digits)
    if (incl_outliers) {
        Q <- round(fivenum(data), digits = digits)
    } else {
        Q <- round(c(boxplot(data, plot = FALSE)$stats), digits = digits)
    }
    res["Min"] <- Q[1]
    res["Q1"]  <- Q[2]
    res["Median"] <- Q[3]
    res["Q3"] <- Q[4]
    res["Max"] <- Q[5]
    return(res)
}

compute_summary.Date <- function(data, group_var = NULL, group = NULL,
                                 incl_outliers, digits) {
    res <- compute_summary.default(unclass(data), group_var, group, incl_outliers,
                                   digits)
    for (i in c("Mean", "Min", "Q1", "Median", "Q3", "Max"))
        class(res[, i]) <- oldClass(data)
    return(res)
}

################################################################################
## Helper for summarize_continuous
prettify.summarize.numeric <- function(x,
                                       colnames = get_option(x, "colnames"),
                                       sep = get_option(x, "sep"),
                                       sanitize = get_option(x, "sanitize"),
                                       count = get_option(x, "count"),
                                       mean_sd = get_option(x, "mean_sd"),
                                       quantiles = get_option(x, "quantiles"),
                                       ...) {

    tab <- x
    ## drop duplicted variable names
    tmp <- tab$variable
    tmp[duplicated(tmp)] <- ""
    tab$variable <- tmp

    ## if not all are TRUE subset results object
    if (!all(count, mean_sd, quantiles)) {

        ## if not any is TRUE (i.e. all are FALSE):
        if (!any(count, mean_sd, quantiles)) {
            stop("Nothing to compute. All quantities are set to FALSE.")
        }
        if (count == FALSE) {
            tab$N <- NULL
            tab$Missing <- NULL
        }
        if (mean_sd == FALSE) {
            tab$Mean <- NULL
            tab$SD <- NULL
        }
        if (quantiles == FALSE) {
            tab$Min <- NULL
            tab$Q1 <- NULL
            tab$Median <- NULL
            tab$Q3 <- NULL
            tab$Max <- NULL
        }
        if (count == FALSE || (mean_sd == FALSE && quantiles == FALSE)) {
            tab$blank_1 <- NULL
        }
        if (mean_sd == FALSE || quantiles == FALSE) {
            tab$blank_2 <- NULL
        }
    }

    if (any(names(tab) == "blank")) {
        start <- which(names(tab) == "blank") + 1
    } else {
        start <- 2
    }
    if (any(grepl("blank_", names(tab)))) {
        idx <- grep("blank_", names(tab))
        if (length(idx) == 1) {
            rules <- paste("  \\cmidrule{", start, "-", idx - 1, "}  ",
                           "\\cmidrule{", idx + 1, "-", length(names(tab)), "}\n",
                           sep = "")
        } else {
            rules <- paste("  \\cmidrule{", start, "-", idx[1] - 1, "}  ",
                           "\\cmidrule{", idx[1] + 1, "-", idx[2] - 1, "} ",
                           "\\cmidrule{", idx[2] + 1, "-", length(names(tab)), "}\n",
                           sep = "")
        }
    } else {
        rules <- paste("  \\cmidrule{", start, "-", length(names(tab)), "}\n",
                       sep = "")
    }

    align <- paste("ll",
                   paste(rep("r", length(names(tab)) - 1), collapse = ""),
                   sep = "")

    ## Define column names
    if (!is.null(colnames)) {
        colNames <- names(tab)
        ## blank doesn't need to be specified in colnames
        if (sum(nms <- !grepl("blank", colNames)) != length(colnames))
            stop(sQuote("colnames"), " has wrong length")
        colNames[nms] <- colnames
    } else {
        colNames <- names(tab)
        colNames[1] <- " "
    }

    colNames[grep("blank", colNames)] <- " "
    colnames(tab) <- colNames

    set_options(tab, colnames = colNames,
                rules = rules, align = align,
                sep = sep, sanitize = sanitize,
                class = "summary")
}


################################################################################
# LaTeX Tables with Descriptves for Factor Variables
summarize_factor <- function(data, variables = names(data),
                             variable.labels = labels, labels = NULL, group = NULL,
                             test = !is.null(group),
                             colnames = NULL, digits = 3, digits.pval = 3,
                             smallest.pval = 0.001, sep = TRUE, sanitize = TRUE,
                             drop = TRUE, show.NAs = any(is.na(data[, variables])),
                             percent = TRUE, cumulative = FALSE,
                             na.lab = "<Missing>", ...) {

    ## get factors
    fac <- mySapply(data[, variables], is.factor)
    ## drop missings
    if (drop) {
        compl.missing <- mySapply(data[, variables], function(x) all(is.na(x)))
        fac <- fac & !compl.missing
    }

    ## if all variables are not factors:
    if (all(!fac))
        stop("None of the variables is a factor or all variables are missing. Nothing to compute.")
    ## if non-factors are dropped:
    if (any(!fac))
        message("Non-factors are dropped from the summary")

    if (is.null(variable.labels)) {
        variable.labels <- variables
    } else {
        if (is.logical(variable.labels) && variable.labels) {
            variable.labels <- labels(data, which = variables)
        } else {
            if (length(variables) != length(variable.labels))
                stop(sQuote("variables"), " and ", sQuote("variable.labels"),
                     " must have the same length")
        }
    }

    ## subset variables to non-factors only
    variables <- variables[fac]
    variable.labels <- variable.labels[fac]

    if (show.NAs) {
        ## convert NAs to factor levels
        if (length(variables) > 1) {
            data[, variables] <- as.data.frame(lapply(data[, variables], NAtoLvl, na.lab = na.lab))
        } else {
            data[, variables] <- NAtoLvl(data[, variables], na.lab)
        }
    }

    if (!is.null(group)) {
        if(!is.factor(data[, group]))
            stop(sQuote("group"), " must be a factor variable")
        if (group %in% variables) {
            idx <- variables != group
            variables <- variables[idx]
            variable.labels <- variable.labels[idx]
        }
        ## remove observations with missing group:
        if (any(is.na(data[, group]))) {
            warning("Removed observations with missing group")
            data <- data[!is.na(data[, group]), ]
        }
        group_var <- data[, group]

        cl <- match.call()
        cl[["group"]] <- NULL
        ## modify call to obtain results for grouped data
        print_single_tabs <- function(level, data, grp_var) {
            dat <- data[grp_var == level, ]
            ## make sure no fatcor level is dropped
            dat <- keep_levels(dat, data)
            cl[["data"]] <- dat
            ## test is not needed in single tables
            cl[["test"]] <- FALSE
            ## set variables
            cl[["variables"]] <- variables
            cl[["prettify"]] <- FALSE
            if (!is.null(variable.labels))
                cl[["variable.labels"]] <- variable.labels
            ## re-evaluate modified call
            eval(cl)
        }
        res <- lapply(levels(group_var), print_single_tabs,
                      data = data, grp_var = group_var)

        res[-1] <- lapply(res[-1], function(x) x[, -c(1:2)])
        stats <- do.call("cbind", res)

        if (!is.character(test) && test)
            test <- "fisher.test"

        if (all(is.character(test))) {
            if (length(test) == 1)
                test <- rep(test, length(variables))
            testdat <- as.matrix(stats[, grep("N", colnames(stats))])
            pval <- rep(NA, length(variables))
            for (i in 1:length(variables)) {
                test_tab <- testdat[stats$variable == unique(stats$variable)[i] & stats$Level != na.lab, ]
                pval[i] <- eval(call(test[i], test_tab))$p.value
            }
            ## make sure rounding is to digits.pval digits
            pval <- format.pval(round(pval, digits = digits.pval),
                                eps = smallest.pval)
            ## make sure not to drop trailing zeros
            pval2 <- suppressWarnings(as.numeric(pval))
            pval[!is.na(pval2)] <- sprintf(paste0("%0.", digits.pval, "f"),
                                           pval2[!is.na(pval2)])
            stats$blank_p <- ""
            stats$p.value <- ""
            stats$p.value[!duplicated(stats$variable)] <- pval
        }

        stats <- set_options(stats, sep = get_option(res[[1]], "sep"),
                             sanitize = get_option(res[[1]], "sanitize"),
                             colnames = get_option(res[[1]], "colnames"),
                             percent = get_option(res[[1]], "percent"),
                             group_labels = paste(group, levels(group_var), sep = ": "),
                             class = "summarize.factor")
        return(prettify(stats))
    }

    ## test not sensible
    if (test || is.character(test))
        warning(sQuote("test"), " is ignored if no ", sQuote("group"), " is given")

    ## repeate variables to match no. of levels:
    n.levels <- mySapply(data[, variables], function(x) length(levels(x)))

    var2 <- unlist(lapply(1:length(variables),
                          function(i) rep(variables[i], each = n.levels[i])))
    var_labels <- unlist(lapply(1:length(variables),
                                function(i) rep(variable.labels[i], each = n.levels[i])))

    ## get all levels
    lvls <- unlist(lapply(variables, function(x) levels(data[, x])))
    colnames(lvls) <- NULL

    ## setup results object
    stats <- data.frame(variable = var_labels, Level = lvls, blank = "",
                        N = NA, Fraction = NA, CumSum = NA,
                        stringsAsFactors = FALSE)
    if (!cumulative) {
        stats$CumSum <- NULL
    }
    rownames(stats) <- NULL

    ## compute statistics
    for (i in 1:length(var2)) {
        notna <- sum(!is.na(data[, var2[i]]))
        stats$N[i] <- sum(data[, var2[i]] == lvls[i], na.rm = TRUE)
        stats$Fraction[i] <- round(stats$N[i]/notna, digits = digits)
        if (cumulative)
            stats$CumSum[i] <- sum(stats$Fraction[1:i][var2[1:i] == var2[i]])
    }
    if (percent) {
        stats$Fraction <- sprintf(paste0("%3.", digits - 2,"f"),
                                  stats$Fraction * 100)
        if (cumulative)
            stats$CumSum <- sprintf(paste0("%3.", digits - 2,"f"),
                                    stats$CumSum * 100)
    } else {
        stats$Fraction <- sprintf(paste0("%1.", digits,"f"), stats$Fraction)
        if (cumulative)
            stats$CumSum <- sprintf(paste0("%1.", digits,"f"), stats$CumSum)
    }
    stats <- set_options(stats, sep = sep, sanitize = sanitize, colnames = colnames,
                         percent = percent, class = "summarize.factor")

    dots <- list(...)
    if (length(dots) > 0 && "prettify" %in% names(dots) && !isTRUE(dots$prettify))
        return(stats)

    prettify(stats)
}

################################################################################
## Helper for summarize_factor
prettify.summarize.factor <- function(x,
                                      colnames = get_option(x, "colnames"),
                                      sep = get_option(x, "sep"),
                                      sanitize = get_option(x, "sanitize"),
                                      ...) {

    tab <- x
    ## drop duplicted variable names
    tmp <- tab$variable
    tmp[duplicated(tmp)] <- ""
    tab$variable <- tmp

    ## define rules
    idx <- c(grep("blank", names(tab)), length(names(tab)) + 1)
    rules <- "  \\cmidrule{2-2} "
    for (i in 1:(length(idx) - 1)) {
        rules <- paste0(rules, "\\cmidrule{", idx[i] + 1, "-", idx[i+1] - 1, "} ")
    }
    rules <- paste0(rules, "\n")

    align <- paste("lll",
                   paste(rep("r", length(names(tab)) - 2), collapse = ""),
                   sep = "")

    ## Define column names
    if (!is.null(colnames)) {
        colNames <- names(tab)
        ## blank doesn't need to be specified in colnames
        if (sum(nms <- !grepl("blank", colNames)) != length(colnames))
            stop(sQuote("colnames"), " has wrong length (should be", sum(nms), ")")
        colNames[nms] <- colnames
    } else {
        colNames <- names(tab)
        if (get_option(x, "percent")) {
            colNames[grepl("Fraction", colNames)] <- "%"
            colNames[grepl("CumSum", colNames)] <- "$\\sum$ %"
        } else {
            colNames[grepl("CumSum", colNames)] <- "$\\sum$"
        }
        colNames[1] <- " "

        header <- ""
        ## if more than one blank add group label
        if (!is.null(get_option(x, "group_labels"))) {
            lab <- get_option(x, "group_labels")
            ## if p.values exist last multicolumn
            ## should not include this column
            if (colNames[length(colNames)] == "p.value")
                idx <- idx[-length(idx)]
            header <- paste(rep("&", idx[1]), collapse = " ")
            for (i in 1:(length(idx) - 1)) {
                header <- paste0(header, "\\multicolumn{", idx[i+1] - idx[i] - 1, "}{c}{", lab[i],"}")
                if (i != length(idx) - 1)
                    header <- paste0(header, " & & ")
            }
            if (colNames[length(colNames)] == "p.value")
                header <- paste0(header, " &  & ")
            header <- paste0(header, "\\\\\n")
        }
    }
    colNames[grep("blank", colNames)] <- " "
    colnames(tab) <- colNames

    set_options(tab, colnames = colNames,
                rules = rules, align = align,
                sep = sep, sanitize = sanitize,
                header = header, group_labels = get_option(x, "group_labels"),
                class = "summary")
}

xtable.summary <- function(x, caption = NULL, label = NULL, align = NULL,
                           digits = NULL, display = NULL, ...) {

    ## options that must be set
    align <- ifelse(is.null(align), get_option(x, "align"), align)

    x <- NextMethod("xtable", object = x, caption = caption, label = label,
                    align = align, digits = digits,
                    display = display, ...)
    class(x) <- c("xtable.summary", "xtable","data.frame")
    return(x)
}

print.xtable.summary <- function(x, rules = NULL, header = NULL,
                                 caption.placement = getOption("xtable.caption.placement", "top"),
                                 hline.after = getOption("xtable.hline.after", NULL),
                                 include.rownames = FALSE,
                                 add.to.row = getOption("xtable.add.to.row", NULL),
                                 booktabs = getOption("xtable.booktabs", TRUE),
                                 sanitize.text.function = get_option(x, "sanitize"),
                                 math.style.negative = getOption("xtable.math.style.negative", TRUE),
                                 math.style.exponents = getOption("xtable.math.style.exponents", TRUE),
                                 tabular.environment = getOption("xtable.tabular.environment", "tabular"),
                                 floating = getOption("xtable.floating", FALSE),
                                 latex.environments = getOption("xtable.latex.environments", c("center")),
                                 ...) {

    ## extract rules and headers from object
    rules <- ifelse(is.null(rules), get_option(x, "rules"), rules)
    tmp <- ifelse(is.null(get_option(x, "header")),
                  "",  get_option(x, "header"))
    header <- ifelse(is.null(header), tmp, header)

    ## add endhead for longtables
    endhead <- ""
    if (tabular.environment == "longtable")
        endhead <- "\\endhead\n"

    if (booktabs)
        cat("%% Output requires \\usepackage{booktabs}.\n")
    if (tabular.environment == "longtable")
        cat("%% Output requires \\usepackage{longtable}.\n")

    ## use centering even if not a float
    if (!floating && latex.environments == "center") {
        cat("\\begin{center}\n")
    }

    ## If caption is given and we don't use a floating environment,
    ## we need to make use of the LaTeX package capt-of
    if (!is.null(caption(x)) && !floating &&
         tabular.environment != "longtable") {
        cat("%% Output requires \\usepackage{capt-of}.\n")
        cat("\\begin{minipage}{\\linewidth}\n",
            "  \\captionof{table}{", caption(x), "}\n",
            ifelse(!is.null(label(x)),
                   paste0("  \\label{", label(x), "}\n"), ""), sep = "")
    }

    ## sanitize object?
    if (is.logical(sanitize.text.function)) {
        if (!sanitize.text.function) {
            sanitize.text.function <- function(x) x
        } else {
            sanitize.text.function <- toLatex
        }
    }

    if (is.null(add.to.row)) {
        if (get_option(x, "sep") == TRUE) {
            pos.rules <- which(x[, 1] != "") - 1
        } else {
            pos.rules <- 0
        }
        ## add endhead to first rule
        midrules <- rep(rules, length(pos.rules))
        midrules[1] <- paste(midrules[1], endhead)
        add.to.row <- list(pos = as.list(c(-1, -1, pos.rules, nrow(x))),
                           command = c("\\toprule\n", header,
                                       midrules, "\\bottomrule\n"))
    }

    if (include.rownames)
        warning(sQuote("include.rownames = TRUE"),
                " is ignored.")

    print.xtable(x,
                 caption.placement = caption.placement,
                 hline.after = hline.after,
                 include.rownames = FALSE,
                 booktabs = booktabs,
                 add.to.row = add.to.row,
                 sanitize.text.function = sanitize.text.function,
                 math.style.negative = math.style.negative,
                 math.style.exponents = math.style.exponents,
                 tabular.environment = tabular.environment,
                 floating = floating,
                 ...)

    if (!is.null(caption(x)) && !floating &&
         tabular.environment != "longtable")
        cat("\\end{minipage}\n")
    ## use centering even if not a float
    if (!floating && latex.environments == "center")
        cat("\\end{center}\n")
    if (!floating && tabular.environment != "longtable"
        && latex.environments != "center")
        cat("\\newline\n")
}


## Add print methods for summarize functions
print.summary <- function(x, ...) {
    if (!is.null(get_option(x, "group_labels"))) {
        cn <- colnames(x)
        lab <- get_option(x, "group_labels")
        x <- rbind(cn, data.frame(x))
        colnames(x) <- rep(" ", ncol(x))
        colnames(x)[cn == "N"] <- lab
        rownames(x) <- c(" ", 1:(nrow(x)-1))
    } else {
        rownames(x) <- NULL
    }
    NextMethod("print")
}
