epitab <-
function (x, y = NULL, method = c("oddsratio", "riskratio", "rateratio"), 
    conf.level = 0.95, rev = c("neither", "rows", "columns", 
        "both"), oddsratio = c("wald", "fisher", "midp", "small"), 
    riskratio = c("wald", "boot", "small"), rateratio = c("wald", 
        "midp"), pvalue = c("fisher.exact", "midp.exact", "chi2"), 
    correction = FALSE, verbose = FALSE) 
{
    method <- match.arg(method)
    if (method == "oddsratio" || method == "riskratio") {
        if (is.matrix(x) && !is.null(y)) {
            stop("y argument should be NULL")
        }
        if (is.null(y)) {
            x <- epitable(x, rev = rev)
        }
        else {
            x <- epitable(x, y, rev = rev)
        }
    }
    if (method == "rateratio") {
        if (is.matrix(x) && !is.null(y)) {
            stop("y argument should be NULL")
        }
        if (is.null(y)) {
            x <- ratetable(x, rev = rev)
        }
        else {
            xn <- substitute(x)
            yn <- substitute(y)
            x <- ratetable(x, y, rev = rev)
            colnames(x) <- c(xn, yn)
        }
    }
    if (method == "oddsratio") {
        oddsratio <- match.arg(oddsratio)
        rr <- oddsratio(x, method = oddsratio, verbose = TRUE, 
            correction = correction, conf.level = conf.level) 
        pvalue <- match.arg(pvalue)
        if (pvalue == "chi2") {
            pval <- rr$p.value[, 3]
        }
        else if (pvalue == "midp.exact") {
            pval <- rr$p.value[, 1]
        }
        else {
            pval <- rr$p.value[, 2]
            pvalue <- "fisher.exact"
        }
        tab <- cbind(rr$x[, 1], rr$p.exp[-nrow(rr$p.exp), 1], 
            rr$x[, 2], rr$p.exp[-nrow(rr$p.exp), 2], rr$measure, 
            pval)
        cn <- colnames(x)
        rownames(tab) <- rownames(x)
        colnames(tab) <- c(cn[1], "p0", cn[2], "p1", "oddsratio", 
            "lower", "upper", "p.value")
        if (!is.null(names(dimnames(x)))) {
            names(dimnames(tab)) <- names(dimnames(x))
        }
        if (verbose) {
            fin <- list(tab = tab, measure = oddsratio, conf.level = conf.level, 
                pvalue = pvalue, x = rr$x, data = rr$data, p.exposed = rr$p.exposed, 
                p.outcome = rr$p.outcome, p.value = rr$p.value, 
                correction = correction)
        }
        else {
            fin <- list(tab = tab, measure = oddsratio, conf.level = conf.level, 
                pvalue = pvalue)
        }
    }
    if (method == "riskratio") {
        riskratio <- match.arg(riskratio)
        rr <- riskratio(x, method = riskratio, verbose = TRUE, 
            correction = correction)
        pvalue <- match.arg(pvalue)
        if (pvalue == "chi2") {
            pval <- rr$p.value[, 3]
        }
        else if (pvalue == "midp.exact") {
            pval <- rr$p.value[, 1]
        }
        else {
            pval <- rr$p.value[, 2]
            pvaue <- "fisher.exact"
        }
        tab <- cbind(rr$x[, 1], rr$p.out[-nrow(rr$p.out), 1], 
            rr$x[, 2], rr$p.out[-nrow(rr$p.out), 2], rr$measure, 
            pval)
        cn <- colnames(x)
        rownames(tab) <- rownames(x)
        colnames(tab) <- c(cn[1], "p0", cn[2], "p1", "riskratio", 
            "lower", "upper", "p.value")
        if (!is.null(names(dimnames(x)))) {
            names(dimnames(tab)) <- names(dimnames(x))
        }
        if (verbose) {
            fin <- list(tab = tab, measure = riskratio, conf.level = conf.level, 
                pvalue = pvalue, x = rr$x, data = rr$data, p.exposed = rr$p.exposed, 
                p.outcome = rr$p.outcome, p.value = rr$p.value, 
                correction = correction)
        }
        else {
            fin <- list(tab = tab, measure = riskratio, conf.level = conf.level, 
                pvalue = pvalue)
        }
    }
    if (method == "rateratio") {
        rateratio <- match.arg(rateratio)
        rr <- rateratio(x, method = rateratio, verbose = TRUE)
        pvalue <- match.arg(pvalue)
        if (pvalue == "chi2") {
            pval <- rr$p.value[, 2]
            pvalue <- "norm.approx"
        }
        else {
            pval <- rr$p.value[, 1]
            pvalue <- "midp.exact"
        }
        tab <- cbind(rr$x, rr$measure, pval)
        cn <- colnames(x)
        rownames(tab) <- rownames(x)
        colnames(tab) <- c(cn[1:2], "rateratio", "lower", "upper", 
            "p.value")
        if (!is.null(names(dimnames(x)))) {
            names(dimnames(tab)) <- names(dimnames(x))
        }
        if (verbose) {
            fin <- list(tab = tab, measure = rateratio, conf.level = conf.level, 
                pvalue = pvalue, x = rr$x, data = rr$data, p.value = rr$p.value)
        }
        else {
            fin <- list(tab = tab, measure = rateratio, conf.level = conf.level, 
                pvalue = pvalue)
        }
    }
    fin
}
