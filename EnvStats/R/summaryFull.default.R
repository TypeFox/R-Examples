summaryFull.default <-
function (object, group = NULL, combine.groups = FALSE, drop.unused.levels = TRUE, 
    rm.group.na = TRUE, stats = NULL, trim = 0.1, sd.method = "sqrt.unbiased", 
    geo.sd.method = "sqrt.unbiased", skew.list = list(), kurtosis.list = list(), 
    cv.list = list(), digits = max(3, getOption("digits") - 3), 
    digit.type = "signif", stats.in.rows = TRUE, drop0trailing = TRUE, 
    data.name = deparse(substitute(object)), ...) 
{
    digit.type <- match.arg(digit.type, c("signif", "round"))
    x <- as.vector(unlist(object))
    if (!is.numeric(x)) 
        stop("All elements of 'object' must be numeric")
    is.null.group <- is.null(group)
    L0 <- !is.null.group && length(group) == 0
    if (L0) 
        warning("length(group) = 0 so the \"group\" argument is ignored.")
    if (!is.null.group & !L0) {
        if (!is.factor(group)) 
            group <- factor(unlist(group))
        n <- length(group)
        if (n != length(x)) 
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
                x <- x[!group.NA]
                warning("Missing values omitted from \"group\" and corresponding elements omitted from \"object\".")
            }
        }
    }
    no.group <- is.null.group || L0 || one.category
    index <- !is.na(x)
    x.le.0 <- sum(index) && any(x[index] <= 0)
    if (is.null(stats)) {
        stats <- ifelse(x.le.0, "for.non.pos", "all")
    }
    else {
        table <- c("all", "for.non.pos", "n", "n.miss", "mean", 
            "median", "trimmed.mean", "geo.mean", "skew", "kurtosis", 
            "min", "max", "range", "1st.quart", "3rd.quart", 
            "sd", "geo.sd", "iqr", "mad", "cv")
        stats <- pmatch(x = stats, table = table, nomatch = 0)
        if (!any(stats > 0) || stats == 0L) 
            stop("invalid value for \"stats\"")
        stats <- table[stats]
        if (identical(stats, "all") && x.le.0) {
            stats <- "for.non.pos"
            warning(paste("'object' contains non-negative values so a geometric mean", 
                "and geometric standard deviation cannot be computed"))
        }
        else if ((!identical(stats, "for.non.pos")) && x.le.0) {
            match.vec <- match(c("geo.mean", "geo.sd"), stats, 
                nomatch = 0)
            if (any(match.vec > 0)) {
                stats <- stats[-match.vec]
                warning(paste("'object' contains non-negative values so a geometric mean", 
                  "and geometric standard deviation cannot be computed"))
            }
        }
    }
    if (no.group) {
        if (length(data.name) > 1) 
            data.name <- "object"
        mat <- summaryFull.vec(x = x, stats = stats, trim = trim, 
            sd.method = sd.method, geo.sd.method = geo.sd.method, 
            skew.list = skew.list, kurtosis.list = kurtosis.list, 
            cv.list = cv.list, digits = digits, digit.type = digit.type, 
            x.name = data.name, stats.in.rows = stats.in.rows)
    }
    else {
        n.groups <- length(levels(group))
        Combined <- summaryFull.vec(x = x, stats = stats, trim = trim, 
            sd.method = sd.method, geo.sd.method = geo.sd.method, 
            skew.list = skew.list, kurtosis.list = kurtosis.list, 
            cv.list = cv.list, digits = digits, digit.type = digit.type, 
            x.name = "Combined", stats.in.rows = TRUE)
        mat <- sapply(split(x, group), summaryFull.vec, stats = stats, 
            trim = trim, sd.method = sd.method, geo.sd.method = geo.sd.method, 
            skew.list = skew.list, kurtosis.list = kurtosis.list, 
            cv.list = cv.list, digits = digits, digit.type = digit.type, 
            stats.in.rows = TRUE)
        rn <- rownames(Combined)
        rownames(mat) <- rn
        if (combine.groups) 
            mat <- cbind(mat, Combined)
        if (!stats.in.rows) 
            mat <- t(mat)
    }
    stat.names.vec <- dimnames(mat)[[ifelse(stats.in.rows, 1, 
        2)]]
    index1 <- match("NA's", stat.names.vec, nomatch = 0)
    index2 <- match("N.Total", stat.names.vec, nomatch = 0)
    if (index1 > 0) {
        if (stats.in.rows) {
            if (all(mat[index1, ] == 0)) 
                mat <- mat[-c(index1, index2), , drop = FALSE]
        }
        else {
            if (all(mat[, index1] == 0)) 
                mat <- mat[, -c(index1, index2), drop = FALSE]
        }
    }
    oldClass(mat) <- "summaryStats"
    attr(mat, "stats.in.rows") <- stats.in.rows
    attr(mat, "drop0trailing") <- as.logical(drop0trailing)
    mat
}
