summaryStats.factor <-
function (object, group = NULL, drop.unused.levels = TRUE, digits = max(3, 
    getOption("digits") - 3), digit.type = "round", drop0trailing = TRUE, 
    show.na = TRUE, show.0.na = FALSE, p.value = FALSE, p.value.digits = 2, 
    p.value.digit.type = "signif", test = "chisq", test.arg.list = NULL, 
    combine.levels = TRUE, combine.groups = FALSE, rm.group.na = TRUE, 
    ci = p.value & test != "chisq", conf.level = 0.95, stats.in.rows = FALSE, 
    ...) 
{
    digit.type <- match.arg(digit.type, c("signif", "round"))
    p.value.digit.type <- match.arg(p.value.digit.type, c("signif", 
        "round"))
    test <- match.arg(test, c("chisq", "prop", "fisher", "binom"))
    if ((p.value | ci) & !combine.levels) {
        combine.levels <- TRUE
        warning("Since p.value=TRUE and/or ci=TRUE, the argument \"combine.levels\" is set to TRUE")
    }
    if (ci & test == "chisq") {
        ci <- FALSE
        warning("test=\"chisq\" so ci set to ci=FALSE")
    }
    is.null.group <- is.null(group)
    L0 <- !is.null.group && length(group) == 0
    if (L0) 
        warning("length(group) = 0 so the \"group\" argument is ignored.")
    if (!is.null.group & !L0) {
        if (!is.factor(group)) 
            group <- factor(unlist(group))
        n <- length(group)
        if (n != length(object)) 
            stop("\"group\" must have the same number of elements as \"object\".")
        group.NA <- is.na(group)
        n.group.NA <- sum(group.NA)
        if (n.group.NA & !rm.group.na) 
            stop("rm.group.na=FALSE so missing values are not allowed in \"group\".")
        if (n.group.NA == n) 
            stop("All values of \"group\" are missing")
        if (n.group.NA) {
            levels.group <- levels(group)
            group <- as.character(group)[!group.NA]
            levels.group <- levels.group[levels.group %in% unique(group)]
            group <- factor(group, levels = levels.group)
        }
        levels.group <- levels(group)
        unique.group <- unique(group)
        if (length(levels.group) > length(unique.group) && drop.unused.levels) 
            group <- factor(group, levels = intersect(levels.group, 
                unique.group))
        one.category <- length(unique.group) == 1
        if (one.category) {
            warning("All non-missing values of \"group\" are the same so this argument is ignored.")
        }
        else {
            if (n.group.NA) {
                levels.object <- levels(object)
                object <- as.character(object)[!group.NA]
                levels.object <- levels.object[levels.object %in% 
                  unique(object)]
                object <- factor(object, levels = levels.object)
                warning("Missing values omitted from \"group\" and corresponding elements omitted from \"object\".")
            }
        }
    }
    levels.object <- levels(object)
    unique.object <- unique(object)
    if (length(levels.object) > length(unique.object) && drop.unused.levels) 
        object <- factor(object, levels = intersect(levels.object, 
            unique.object))
    n.levels <- length(levels(object))
    if (n.levels < 2 & (p.value | ci)) {
        p.value <- FALSE
        ci <- FALSE
        warning("Number of categories in \"object\" less than 2 so p-value and/or CI cannot be computed")
    }
    if (is.null.group || L0 || one.category) {
        mat <- summaryStats.factor.vec(x = object, combine.levels = combine.levels, 
            digits = digits, digit.type = digit.type, show.na = show.na, 
            show.0.na = show.0.na)
        if (p.value || ci) {
            if (test == "fisher") 
                stop("test=\"fisher\" not available when there is no \"group\" argument")
            if (n.levels > 2) {
                if (test == "prop") 
                  stop("test=\"prop\" only available with 2 categories")
                if (test == "binom") 
                  stop("test=\"binom\" only available with 2 categories")
            }
            tab.object <- table(object)
            switch(test, chisq = {
                test.list <- do.call("chisq.test", args = c(list(x = tab.object), 
                  test.arg.list))
                p.col.name <- "ChiSq_p"
            }, prop = {
                test.list <- do.call("prop.test", args = c(list(x = tab.object[1], 
                  n = sum(tab.object)), test.arg.list))
                p.col.name <- "ChiSq_p"
            }, binom = {
                test.list <- do.call("binom.test", args = c(list(x = tab.object[1], 
                  n = sum(tab.object)), test.arg.list))
                p.col.name <- "Exact_p"
            })
            if (p.value) {
                p <- test.list$p.value
                p.column <- rep(NA, nrow(mat))
                rn <- rownames(mat)
                cn <- colnames(mat)
                names(p.column) <- rn
                p.column[match("Combined", rn)] <- do.call(p.value.digit.type, 
                  args = list(x = p, digits = p.value.digits))
                mat <- cbind(mat, p.column)
                colnames(mat) <- c(cn, p.col.name)
            }
            if (ci & n.levels == 2) {
                conf.int <- 100 * test.list$conf.int
                LCL.column <- rep(NA, nrow(mat))
                LCL.column[1] <- do.call(digit.type, args = list(x = conf.int[1], 
                  digits = digits))
                UCL.column <- rep(NA, nrow(mat))
                UCL.column[1] <- do.call(digit.type, args = list(x = conf.int[2], 
                  digits = digits))
                rn <- rownames(mat)
                cn <- colnames(mat)
                names(LCL.column) <- rn
                names(UCL.column) <- rn
                mat <- cbind(mat, LCL.column, UCL.column)
                colnames(mat) <- c(cn, paste(100 * conf.level, 
                  "%.", c("LCL", "UCL"), sep = ""))
            }
        }
    }
    else {
        if (p.value & test == "binom") 
            stop("You cannot set test=\"binom\" when \"group\" is supplied")
        levels.group <- levels(group)
        dum.list <- lapply(split(object, group), summaryStats.factor.vec, 
            combine.levels = combine.levels, digits = digits, 
            digit.type = digit.type, show.na = show.na, show.0.na = TRUE)
        mat <- do.call("cbind", dum.list)
        colnames(mat) <- paste(rep(levels.group, each = 2), "(", 
            colnames(mat), ")", sep = "")
        if (show.na && !show.0.na) {
            rn <- rownames(mat)
            index <- match("NA's", rn)
            na.row <- mat[index, ]
            if (all(na.row[!is.na(na.row)] == 0)) 
                mat <- mat[-index, ]
        }
        if (combine.groups) {
            Combined <- summaryStats.factor.vec(object, combine.levels = combine.levels, 
                digits = digits, digit.type = digit.type, show.na = show.na, 
                show.0.na = show.0.na)
            colnames(Combined) <- paste("Combined", "(", colnames(Combined), 
                ")", sep = "")
            mat <- cbind(mat, Combined)
        }
        if (p.value | ci) {
            test.tab <- table(object, group)
            dim.test.tab <- dim(test.tab)
            if (any(dim.test.tab > 2) & test == "prop") 
                stop("When \"object\" and/or \"group\" has more than 2 levels you cannot set test=\"prop\"")
            if (test == "prop") 
                test.list <- do.call("prop.test", args = c(list(x = test.tab[1, 
                  ], n = colSums(test.tab)), test.arg.list))
            else test.list <- do.call(paste(test, "test", sep = "."), 
                args = c(list(x = test.tab), test.arg.list))
            if (p.value) {
                p <- test.list$p.value
                col.name <- ifelse(test == "chisq" | test == 
                  "prop", "ChiSq_p", "Fisher_p")
                p.column <- rep(NA, nrow(mat))
                rn <- rownames(mat)
                cn <- colnames(mat)
                names(p.column) <- rn
                p.column[match("Combined", rn)] <- do.call(p.value.digit.type, 
                  args = list(x = p, digits = p.value.digits))
                mat <- cbind(mat, p.column)
                colnames(mat) <- c(cn, col.name)
            }
            if (ci & all(dim.test.tab == 2)) {
                conf.int <- test.list$conf.int
                if (test != "fisher") 
                  conf.int <- 100 * conf.int
                col.names <- paste(100 * conf.level, "%.", c("LCL", 
                  "UCL"), sep = "")
                if (test == "fisher") 
                  col.names <- paste(col.names, "OR", sep = ".")
                else if (test == "prop") 
                  col.names <- paste(col.names, "between", sep = ".")
                LCL.column <- rep(NA, nrow(mat))
                UCL.column <- rep(NA, nrow(mat))
                rn <- rownames(mat)
                cn <- colnames(mat)
                names(LCL.column) <- rn
                names(UCL.column) <- rn
                LCL.column[match("Combined", rn)] <- do.call(digit.type, 
                  args = list(x = conf.int[1], digits = digits))
                UCL.column[match("Combined", rn)] <- do.call(digit.type, 
                  args = list(x = conf.int[2], digits = digits))
                mat <- cbind(mat, LCL.column, UCL.column)
                colnames(mat) <- c(cn, col.names)
            }
        }
    }
    if (stats.in.rows) 
        mat <- t(mat)
    oldClass(mat) <- "summaryStats"
    attr(mat, "stats.in.rows") <- stats.in.rows
    attr(mat, "drop0trailing") <- as.logical(drop0trailing)
    mat
}
