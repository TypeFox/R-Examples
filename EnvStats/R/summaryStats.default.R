summaryStats.default <-
function (object, group = NULL, drop.unused.levels = TRUE, se = FALSE, 
    quartiles = FALSE, digits = max(3, getOption("digits") - 
        3), digit.type = "round", drop0trailing = TRUE, show.na = TRUE, 
    show.0.na = FALSE, p.value = FALSE, p.value.digits = 2, p.value.digit.type = "signif", 
    test = "parametric", test.arg.list = NULL, combine.groups = p.value, 
    rm.group.na = TRUE, group.p.value.type = NULL, alternative = "two.sided", 
    ci = NULL, ci.between = NULL, conf.level = 0.95, stats.in.rows = FALSE, 
    data.name = deparse(substitute(object)), ...) 
{
    digit.type <- match.arg(digit.type, c("signif", "round"))
    p.value.digit.type <- match.arg(p.value.digit.type, c("signif", 
        "round"))
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
    if (is.null(group.p.value.type)) 
        group.p.value.type <- ifelse(combine.groups, "between", 
            "within")
    group.p.value.type <- match.arg(group.p.value.type, c("between", 
        "within"))
    if (is.null(ci)) 
        ci <- p.value & (no.group | (!no.group & group.p.value.type == 
            "within"))
    if (is.null(ci.between)) 
        ci.between <- p.value & group.p.value.type == "between"
    test <- match.arg(test, c("parametric", "nonparametric"))
    alternative <- match.arg(alternative, c("two.sided", "less", 
        "greater"))
    if (no.group) {
        if (length(data.name) > 1) 
            data.name <- "object"
        mat <- summaryStats.vec(x = x, digits = digits, digit.type = digit.type, 
            se = se, quartiles = quartiles, show.na = show.na, 
            show.0.na = show.0.na, p.value = p.value, p.value.digits = p.value.digits, 
            p.value.digit.type = p.value.digit.type, test = test, 
            test.arg.list = test.arg.list, alternative = alternative, 
            ci = ci, conf.level = conf.level, x.name = data.name, 
            stats.in.rows = stats.in.rows)
    }
    else {
        n.groups <- length(levels(group))
        if (p.value && group.p.value.type == "between" && !combine.groups) {
            combine.groups <- TRUE
            warning("group.p.value.type=\"between\" so \"combine.groups\" set to TRUE")
        }
        Combined <- summaryStats.vec(x = x, digits = digits, 
            digit.type = digit.type, se = se, quartiles = quartiles, 
            show.na = show.na, show.0.na = TRUE, p.value = p.value & 
                group.p.value.type == "within", p.value.digits = p.value.digits, 
            p.value.digit.type = p.value.digit.type, test = test, 
            test.arg.list = test.arg.list, alternative = alternative, 
            ci = ci, conf.level = conf.level, x.name = "Combined", 
            stats.in.rows = TRUE)
        mat <- sapply(split(x, group), summaryStats.vec, digits = digits, 
            digit.type = digit.type, se = se, quartiles = quartiles, 
            show.na = show.na, show.0.na = TRUE, p.value = p.value & 
                group.p.value.type == "within", p.value.digits = p.value.digits, 
            p.value.digit.type = p.value.digit.type, test = test, 
            test.arg.list = test.arg.list, alternative = alternative, 
            ci = ci, conf.level = conf.level, stats.in.rows = TRUE)
        rn <- rownames(Combined)
        rownames(mat) <- rn
        if (combine.groups) {
            SDs <- mat["SD", ]
            mat <- cbind(mat, Combined)
            if (p.value && group.p.value.type == "between") {
                cn <- colnames(mat)
                nr <- nrow(mat)
                nc <- ncol(mat)
                if (n.groups == 2) {
                  x.1 <- x[group == levels(group)[1]]
                  x.2 <- x[group == levels(group)[2]]
                  if (any(SDs > 0)) {
                    if (test == "parametric") {
                      if (is.null(test.arg.list) || all(is.na(pmatch(names(test.arg.list), 
                        "var.equal")))) 
                        test.arg.list <- c(test.arg.list, list(var.equal = TRUE))
                      test.list <- do.call("t.test", args = c(list(x = x.2, 
                        y = x.1, alternative = alternative, conf.level = conf.level), 
                        test.arg.list))
                      diff.locations <- -diff(test.list$estimate)
                    }
                    else {
                      test.list <- do.call("wilcox.test", args = c(list(x = x.2, 
                        y = x.1, alternative = alternative, conf.int = TRUE, 
                        conf.level = conf.level), test.arg.list))
                      diff.locations <- test.list$estimate
                    }
                    p <- test.list$p.value
                    conf.int <- test.list$conf.int
                  }
                  else {
                    diff.locations <- NA
                    p <- NA
                    conf.int <- c(NA, NA)
                  }
                  LCL <- c(rep(NA, nc - 1), do.call(digit.type, 
                    list(x = conf.int[1], digits = digits)))
                  names(LCL) <- cn
                  UCL <- c(rep(NA, nc - 1), do.call(digit.type, 
                    list(x = conf.int[2], digits = digits)))
                  names(UCL) <- cn
                }
                else {
                  if (any(SDs > 0)) {
                    if (test == "parametric") 
                      p <- anova(lm(x ~ group))[1, "Pr(>F)"]
                    else p <- kruskal.test(x ~ group)$p.value
                  }
                  else p <- NA
                }
                if (n.groups == 2) {
                  diff.locations <- c(rep(NA, nc - 1), do.call(digit.type, 
                    list(x = diff.locations, digits = digits)))
                  names(diff.locations) <- cn
                  index <- match("NA's", rn)
                  if (!is.na(index)) {
                    mat <- rbind(mat[1:(index - 1), ], Diff = diff.locations, 
                      mat[index:nr, ])
                  }
                  else mat <- rbind(mat, Diff = diff.locations)
                  nr <- nrow(mat)
                  rn <- rownames(mat)
                }
                p <- c(rep(NA, nc - 1), do.call(p.value.digit.type, 
                  list(x = p, digits = p.value.digits)))
                names(p) <- cn
                index <- match("NA's", rn)
                if (!is.na(index)) {
                  mat <- rbind(mat[1:(index - 1), ], p.value.between = p, 
                    mat[index:nr, ])
                }
                else mat <- rbind(mat, p.value.between = p)
                nr <- nrow(mat)
                if (n.groups == 2) {
                  if (ci.between) {
                    index <- match("p.value.between", rownames(mat))
                    if (index == nr) 
                      mat <- rbind(mat, LCL, UCL)
                    else mat <- rbind(mat[1:index, ], LCL, UCL, 
                      mat[(index + 1):nr, ])
                    rownames(mat)[index + (1:2)] <- paste(100 * 
                      conf.level, "%.", c("LCL", "UCL"), ".between", 
                      sep = "")
                  }
                }
                string <- "p.value.between"
                if (test == "nonparametric") {
                  string2 <- ifelse(n.groups == 2, "Wilcoxon", 
                    "Kruskal")
                  string <- paste(string2, string, sep = ".")
                }
                else {
                  if (n.groups == 2) {
                    if (is.null(test.arg.list) || all(is.na(pmatch(names(test.arg.list), 
                      "var.equal")))) 
                      var.equal <- TRUE
                    else var.equal <- unlist(test.arg.list[!is.na(pmatch(names(test.arg.list), 
                      "var.equal"))])
                    if (!var.equal) 
                      string <- paste("Welch", string, sep = ".")
                  }
                }
                if (alternative != "two.sided") 
                  string <- paste(string, alternative, sep = ".")
                index <- match("p.value.between", rownames(mat))
                rownames(mat)[index] <- string
            }
        }
        if (show.na && !show.0.na) {
            rn <- rownames(mat)
            index <- match("NA's", rn)
            if (all(mat[index, ] == 0)) 
                mat <- mat[-(index:(index + 1)), ]
        }
        if (!stats.in.rows) 
            mat <- t(mat)
    }
    oldClass(mat) <- "summaryStats"
    attr(mat, "stats.in.rows") <- stats.in.rows
    attr(mat, "drop0trailing") <- as.logical(drop0trailing)
    mat
}
