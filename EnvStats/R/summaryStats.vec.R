summaryStats.vec <-
function (x, digits = max(3, getOption("digits") - 3), digit.type = "round", 
    se = FALSE, quartiles = FALSE, show.na = FALSE, show.0.na = FALSE, 
    p.value = FALSE, p.value.digits = 2, p.value.digit.type = "signif", 
    test = "parametric", test.arg.list = NULL, alternative = "two.sided", 
    ci = FALSE, conf.level = 0.95, x.name = NULL, stats.in.rows = TRUE) 
{
    digit.type <- match.arg(digit.type, c("signif", "round"))
    p.value.digit.type <- match.arg(p.value.digit.type, c("signif", 
        "round"))
    if (is.null(x.name)) 
        x.name <- deparse(substitute(x))
    x <- as.vector(unlist(x))
    if (length(x) < 1 || !is.numeric(x)) {
        x <- numeric(0)
        warning("length(x) < 1 and/or x not numeric.")
    }
    is.na.x <- is.na(x)
    n <- sum(!is.na.x)
    n.na <- sum(is.na.x)
    x <- x[!is.na.x]
    if (n > 0) {
        sd.x <- stats::sd(x)
        vec <- c(N = n, Mean = mean(x), SD = sd.x, Median = median(x), 
            Min = min(x), Max = max(x))
        if (se) {
            vec <- c(vec[c("N", "Mean", "SD")], SE = sd.x/sqrt(n), 
                vec[c("Median", "Min", "Max")])
        }
        if (quartiles) {
            quarts <- stats::quantile(x, probs = c(0.25, 0.75))
            names(quarts) <- c("1st Qu.", "3rd Qu.")
            vec <- c(vec, quarts)
        }
        vec <- do.call(digit.type, list(x = vec, digits = digits))
        if (p.value || ci) {
            if (n == 1 || !(sd.x > 0)) {
                if (p.value) {
                  vec <- c(vec, p.value = NA)
                }
                if (ci) {
                  conf.int <- c(NA, NA)
                  names(conf.int) <- paste(100 * conf.level, 
                    "%.", c("LCL", "UCL"), sep = "")
                  vec <- c(vec, conf.int)
                }
            }
            else {
                if (test == "parametric") {
                  dum.list <- do.call("t.test", args = c(list(x = x, 
                    alternative = alternative), test.arg.list))
                }
                else {
                  dum.list <- do.call("wilcox.test", args = c(list(x = x, 
                    alternative = alternative, conf.int = ci), 
                    test.arg.list))
                }
                p.val <- dum.list$p.value
                if (p.value) {
                  vec <- c(vec, p.value = do.call(p.value.digit.type, 
                    list(x = p.val, digits = p.value.digits)))
                  if (test == "nonparametric") {
                    names(vec)[length(names(vec))] <- "Wilcoxon.p.value"
                  }
                }
                if (ci) {
                  conf.int <- dum.list$conf.int
                  names(conf.int) <- paste(100 * conf.level, 
                    "%.", c("LCL", "UCL"), sep = "")
                  if (test == "nonparametric") {
                    names(conf.int) <- paste("Wilcoxon", names(conf.int), 
                      sep = ".")
                  }
                  vec <- c(vec, do.call(digit.type, list(x = conf.int, 
                    digits = digits)))
                }
            }
        }
    }
    else {
        if (se) 
            vec <- c(N = n, Mean = NA, SD = NA, SE = NA, Median = NA, 
                Min = NA, Max = NA)
        else vec <- c(N = n, Mean = NA, SD = NA, Median = NA, 
            Min = NA, Max = NA)
        if (quartiles) 
            vec <- c(vec, `1st Qu.` = NA, `3rd Qu.` = NA)
        if (p.value) 
            vec <- c(vec, p.value = NA)
        if (ci) {
            conf.int <- c(LCL = NA, UCL = NA)
            names(conf.int) <- paste(100 * conf.level, "%.", 
                c("LCL", "UCL"), sep = "")
            vec <- c(vec, conf.int)
        }
    }
    if (show.na && (n.na > 0 || show.0.na)) {
        vec <- c(vec, `NA's` = n.na, N.Total = n + n.na)
    }
    ret.val <- t(vec)
    dimnames(ret.val)[[1]] <- x.name
    if (stats.in.rows) {
        ret.val <- t(ret.val)
    }
    ret.val
}
