rma.peto <-
function (ai, bi, ci, di, n1i, n2i, data, slab, subset, add = 1/2, 
    to = "only0", drop00 = TRUE, level = 95, digits = 4, verbose = FALSE) 
{
    if (length(add) == 1) 
        add <- c(add, 0)
    if (length(add) != 2) 
        stop("Argument 'add' should specify one or two values (see 'help(rma.peto)').")
    if (length(to) == 1) 
        to <- c(to, "none")
    if (length(to) != 2) 
        stop("Argument 'to' should specify one or two values (see 'help(rma.peto)').")
    if (length(drop00) == 1) 
        drop00 <- c(drop00, FALSE)
    if (length(drop00) != 2) 
        stop("Argument 'drop00' should specify one or two values (see 'help(rma.peto)').")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    if (!is.element(to[1], c("all", "only0", "if0all", "none"))) 
        stop("Unknown 'to' argument specified.")
    if (!is.element(to[2], c("all", "only0", "if0all", "none"))) 
        stop("Unknown 'to' argument specified.")
    measure <- "PETO"
    if (verbose) 
        message("Extracting data and computing yi/vi values ...")
    if (missing(data)) 
        data <- NULL
    if (is.null(data)) {
        data <- sys.frame(sys.parent())
    }
    else {
        if (!is.data.frame(data)) {
            data <- data.frame(data)
        }
    }
    mf <- match.call()
    mf.slab <- mf[[match("slab", names(mf))]]
    mf.subset <- mf[[match("subset", names(mf))]]
    slab <- eval(mf.slab, data, enclos = sys.frame(sys.parent()))
    subset <- eval(mf.subset, data, enclos = sys.frame(sys.parent()))
    mf.ai <- mf[[match("ai", names(mf))]]
    mf.bi <- mf[[match("bi", names(mf))]]
    mf.ci <- mf[[match("ci", names(mf))]]
    mf.di <- mf[[match("di", names(mf))]]
    mf.n1i <- mf[[match("n1i", names(mf))]]
    mf.n2i <- mf[[match("n2i", names(mf))]]
    ai <- eval(mf.ai, data, enclos = sys.frame(sys.parent()))
    bi <- eval(mf.bi, data, enclos = sys.frame(sys.parent()))
    ci <- eval(mf.ci, data, enclos = sys.frame(sys.parent()))
    di <- eval(mf.di, data, enclos = sys.frame(sys.parent()))
    n1i <- eval(mf.n1i, data, enclos = sys.frame(sys.parent()))
    n2i <- eval(mf.n2i, data, enclos = sys.frame(sys.parent()))
    if (is.null(bi)) 
        bi <- n1i - ai
    if (is.null(di)) 
        di <- n2i - ci
    ni <- ai + bi + ci + di
    k <- length(ai)
    ids <- seq_len(k)
    if (verbose) 
        message("Generating/extracting study labels ...")
    if (is.null(slab)) {
        slab.null <- TRUE
        slab <- ids
    }
    else {
        if (anyNA(slab)) 
            stop("NAs in study labels.")
        if (length(slab) != k) 
            stop("Study labels not of same length as data.")
        slab.null <- FALSE
    }
    if (!is.null(subset)) {
        if (verbose) 
            message("Subsetting ...")
        ai <- ai[subset]
        bi <- bi[subset]
        ci <- ci[subset]
        di <- di[subset]
        ni <- ni[subset]
        slab <- slab[subset]
        ids <- ids[subset]
        k <- length(ai)
    }
    if (anyDuplicated(slab)) 
        slab <- make.unique(as.character(slab))
    dat <- escalc(measure = "PETO", ai = ai, bi = bi, ci = ci, 
        di = di, add = add[1], to = to[1], drop00 = drop00[1])
    yi <- dat$yi
    vi <- dat$vi
    if (drop00[2]) {
        id00 <- c(ai == 0L & ci == 0L) | c(bi == 0L & di == 0L)
        id00[is.na(id00)] <- FALSE
        ai[id00] <- NA
        bi[id00] <- NA
        ci[id00] <- NA
        di[id00] <- NA
    }
    ai.f <- ai
    bi.f <- bi
    ci.f <- ci
    di.f <- di
    yi.f <- yi
    vi.f <- vi
    ni.f <- ni
    k.f <- k
    aibicidi.na <- is.na(ai) | is.na(bi) | is.na(ci) | is.na(di)
    if (any(aibicidi.na)) {
        if (verbose) 
            message("Handling NAs in table data ...")
        not.na <- !aibicidi.na
        if (na.act == "na.omit" || na.act == "na.exclude" || 
            na.act == "na.pass") {
            ai <- ai[not.na]
            bi <- bi[not.na]
            ci <- ci[not.na]
            di <- di[not.na]
            k <- length(ai)
            warning("Tables with NAs omitted from model fitting.")
        }
        if (na.act == "na.fail") 
            stop("Missing values in tables.")
    }
    else {
        not.na <- rep(TRUE, k)
    }
    if (k < 1) 
        stop("Processing terminated since k = 0.")
    yivi.na <- is.na(yi) | is.na(vi)
    if (any(yivi.na)) {
        if (verbose) 
            message("Handling NAs in yi/vi ...")
        not.na.yivi <- !yivi.na
        if (na.act == "na.omit" || na.act == "na.exclude" || 
            na.act == "na.pass") {
            yi <- yi[not.na.yivi]
            vi <- vi[not.na.yivi]
            ni <- ni[not.na.yivi]
            warning("Some yi/vi values are NA.")
            attr(yi, "measure") <- measure
            attr(yi, "ni") <- ni
        }
        if (na.act == "na.fail") 
            stop("Missing yi/vi values.")
    }
    else {
        not.na.yivi <- rep(TRUE, k)
    }
    k.yi <- length(yi)
    if (to[2] == "all") {
        ai <- ai + add[2]
        bi <- bi + add[2]
        ci <- ci + add[2]
        di <- di + add[2]
    }
    if (to[2] == "only0") {
        id0 <- c(ai == 0L | bi == 0L | ci == 0L | di == 0L)
        ai[id0] <- ai[id0] + add[2]
        bi[id0] <- bi[id0] + add[2]
        ci[id0] <- ci[id0] + add[2]
        di[id0] <- di[id0] + add[2]
    }
    if (to[2] == "if0all") {
        id0 <- c(ai == 0L | bi == 0L | ci == 0L | di == 0L)
        if (any(id0)) {
            ai <- ai + add[2]
            bi <- bi + add[2]
            ci <- ci + add[2]
            di <- di + add[2]
        }
    }
    n1i <- ai + bi
    n2i <- ci + di
    Ni <- ai + bi + ci + di
    alpha <- ifelse(level > 1, (100 - level)/100, 1 - level)
    if (verbose) 
        message("Model fitting ...")
    xt <- ai + ci
    yt <- bi + di
    Ei <- xt * n1i/Ni
    Vi <- xt * yt * (n1i/Ni) * (n2i/Ni)/(Ni - 1)
    sumVi <- sum(Vi)
    if (sumVi == 0L) 
        stop("One of the two outcomes never occurred in any of the tables. Peto's method cannot be used.")
    b <- sum(ai - Ei)/sumVi
    se <- sqrt(1/sumVi)
    zval <- b/se
    pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
    ci.lb <- b - qnorm(alpha/2, lower.tail = FALSE) * se
    ci.ub <- b + qnorm(alpha/2, lower.tail = FALSE) * se
    names(b) <- "intrcpt"
    vb <- matrix(se^2, dimnames = list("intrcpt", "intrcpt"))
    if (verbose) 
        message("Heterogeneity testing ...")
    k.pos <- sum(Vi > 0)
    Vi[Vi == 0] <- NA
    QE <- sum((ai - Ei)^2/Vi, na.rm = TRUE) - sum(ai - Ei)^2/sum(Vi, 
        na.rm = TRUE)
    if (k.pos >= 1) {
        QEp <- pchisq(QE, df = k.yi - 1, lower.tail = FALSE)
    }
    else {
        QEp <- 1
    }
    wi <- 1/vi
    RSS <- sum(wi * (yi - b)^2)
    if (verbose) 
        message("Computing fit statistics and log likelihood ...")
    ll.ML <- -1/2 * (k.yi) * log(2 * base::pi) - 1/2 * sum(log(vi)) - 
        1/2 * RSS
    ll.REML <- -1/2 * (k.yi - 1) * log(2 * base::pi) + 1/2 * 
        log(k.yi) - 1/2 * sum(log(vi)) - 1/2 * log(sum(wi)) - 
        1/2 * RSS
    dev.ML <- -2 * (ll.ML - sum(dnorm(yi, mean = yi, sd = sqrt(vi), 
        log = TRUE)))
    AIC.ML <- -2 * ll.ML + 2
    BIC.ML <- -2 * ll.ML + log(k.yi)
    AICc.ML <- -2 * ll.ML + 2 * max(k.yi, 3)/(max(k.yi, 3) - 
        2)
    dev.REML <- -2 * (ll.REML - 0)
    AIC.REML <- -2 * ll.REML + 2
    BIC.REML <- -2 * ll.REML + log(k.yi - 1)
    AICc.REML <- -2 * ll.REML + 2 * max(k.yi - 1, 3)/(max(k.yi - 
        1, 3) - 2)
    fit.stats <- matrix(c(ll.ML, dev.ML, AIC.ML, BIC.ML, AICc.ML, 
        ll.REML, dev.REML, AIC.REML, BIC.REML, AICc.REML), ncol = 2, 
        byrow = FALSE)
    dimnames(fit.stats) <- list(c("ll", "dev", "AIC", "BIC", 
        "AICc"), c("ML", "REML"))
    fit.stats <- data.frame(fit.stats)
    if (verbose) 
        message("Preparing output ...")
    parms <- 1
    p <- 1
    p.eff <- 1
    k.eff <- k
    tau2 <- 0
    X.f <- cbind(rep(1, k.f))
    intercept <- TRUE
    int.only <- TRUE
    method <- "FE"
    weighted <- TRUE
    knha <- FALSE
    dfs <- NA
    res <- list(b = b, se = se, zval = zval, pval = pval, ci.lb = ci.lb, 
        ci.ub = ci.ub, vb = vb, tau2 = tau2, k = k, k.f = k.f, 
        k.yi = k.yi, k.pos = k.pos, k.eff = k.eff, p = p, parms = parms, 
        QE = QE, QEp = QEp, int.only = int.only, yi = yi, vi = vi, 
        yi.f = yi.f, vi.f = vi.f, X.f = X.f, ai = ai, bi = bi, 
        ci = ci, di = di, ai.f = ai.f, bi.f = bi.f, ci.f = ci.f, 
        di.f = di.f, ni = ni, ni.f = ni.f, ids = ids, not.na = not.na, 
        not.na.yivi = not.na.yivi, slab = slab, slab.null = slab.null, 
        measure = measure, method = method, weighted = weighted, 
        knha = knha, dfs = dfs, intercept = intercept, digits = digits, 
        level = level, add = add, to = to, drop00 = drop00, fit.stats = fit.stats, 
        call = mf)
    class(res) <- c("rma.peto", "rma")
    return(res)
}
