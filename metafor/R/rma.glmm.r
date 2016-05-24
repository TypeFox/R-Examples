rma.glmm <-
function (ai, bi, ci, di, n1i, n2i, x1i, x2i, t1i, t2i, xi, mi, 
    ti, ni, mods, measure, intercept = TRUE, data, slab, subset, 
    add = 1/2, to = "only0", drop00 = TRUE, vtype = "LS", model = "UM.FS", 
    method = "ML", tdist = FALSE, level = 95, digits = 4, btt, 
    nAGQ = 7, verbose = FALSE, control) 
{
    if (missing(measure)) 
        stop("Need to specify 'measure' argument.")
    if (!is.element(measure, c("OR", "IRR", "PLO", "IRLN"))) 
        stop("Unknown 'measure' specified.")
    if (!is.element(method, c("FE", "ML"))) 
        stop("Unknown 'method' specified.")
    if (length(add) > 1) 
        add <- add[1]
    if (length(to) > 1) 
        to <- to[1]
    if (!is.element(model, c("UM.FS", "UM.RS", "CM.EL", "CM.AL"))) 
        stop("Unknown 'model' specified.")
    if (model == "CM.AL" && measure == "IRR") 
        model <- "CM.EL"
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    if (missing(btt)) 
        btt <- NULL
    if (missing(control)) 
        control <- list()
    if (is.element(measure, c("OR", "IRR")) && model == "UM.RS" && 
        method == "ML" && nAGQ > 1) {
        warning("Currently not possible to fit RE/ME model='UM.RS' with nAGQ > 1. nAGQ automatically set to 1.")
        nAGQ <- 1
    }
    knha <- tdist
    very.verbose <- ifelse(!is.logical(verbose) && verbose > 
        1, TRUE, FALSE)
    if (very.verbose) 
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
    mf.mods <- mf[[match("mods", names(mf))]]
    slab <- eval(mf.slab, data, enclos = sys.frame(sys.parent()))
    subset <- eval(mf.subset, data, enclos = sys.frame(sys.parent()))
    mods <- eval(mf.mods, data, enclos = sys.frame(sys.parent()))
    ai <- bi <- ci <- di <- x1i <- x2i <- t1i <- t2i <- xi <- mi <- ti <- ni <- NA
    if (is.element(measure, c("OR"))) {
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
        k <- length(ai)
        if (!is.null(subset)) {
            ai <- ai[subset]
            bi <- bi[subset]
            ci <- ci[subset]
            di <- di[subset]
        }
        dat <- escalc(measure = measure, ai = ai, bi = bi, ci = ci, 
            di = di, add = add, to = to, drop00 = drop00, vtype = vtype)
    }
    if (is.element(measure, c("IRR"))) {
        mf.x1i <- mf[[match("x1i", names(mf))]]
        mf.x2i <- mf[[match("x2i", names(mf))]]
        mf.t1i <- mf[[match("t1i", names(mf))]]
        mf.t2i <- mf[[match("t2i", names(mf))]]
        x1i <- eval(mf.x1i, data, enclos = sys.frame(sys.parent()))
        x2i <- eval(mf.x2i, data, enclos = sys.frame(sys.parent()))
        t1i <- eval(mf.t1i, data, enclos = sys.frame(sys.parent()))
        t2i <- eval(mf.t2i, data, enclos = sys.frame(sys.parent()))
        k <- length(x1i)
        if (!is.null(subset)) {
            x1i <- x1i[subset]
            x2i <- x2i[subset]
            t1i <- t1i[subset]
            t2i <- t2i[subset]
        }
        dat <- escalc(measure = measure, x1i = x1i, x2i = x2i, 
            t1i = t1i, t2i = t2i, add = add, to = to, drop00 = drop00, 
            vtype = vtype)
    }
    if (is.element(measure, c("PLO"))) {
        mf.xi <- mf[[match("xi", names(mf))]]
        mf.mi <- mf[[match("mi", names(mf))]]
        mf.ni <- mf[[match("ni", names(mf))]]
        xi <- eval(mf.xi, data, enclos = sys.frame(sys.parent()))
        mi <- eval(mf.mi, data, enclos = sys.frame(sys.parent()))
        ni <- eval(mf.ni, data, enclos = sys.frame(sys.parent()))
        if (is.null(mi)) 
            mi <- ni - xi
        k <- length(xi)
        if (!is.null(subset)) {
            xi <- xi[subset]
            mi <- mi[subset]
        }
        dat <- escalc(measure = measure, xi = xi, mi = mi, add = add, 
            to = to, vtype = vtype)
    }
    if (is.element(measure, c("IRLN"))) {
        mf.xi <- mf[[match("xi", names(mf))]]
        mf.ti <- mf[[match("ti", names(mf))]]
        xi <- eval(mf.xi, data, enclos = sys.frame(sys.parent()))
        ti <- eval(mf.ti, data, enclos = sys.frame(sys.parent()))
        k <- length(xi)
        if (!is.null(subset)) {
            xi <- xi[subset]
            ti <- ti[subset]
        }
        dat <- escalc(measure = measure, xi = xi, ti = ti, add = add, 
            to = to, vtype = vtype)
    }
    yi <- dat$yi
    vi <- dat$vi
    ni <- attr(yi, "ni")
    ids <- seq_len(k)
    if (very.verbose) 
        message("Creating model matrix ...")
    is.formula <- class(mods) == "formula"
    if (is.formula) {
        options(na.action = "na.pass")
        mods <- model.matrix(mods, data = data)
        attr(mods, "assign") <- NULL
        options(na.action = na.act)
        intercept <- FALSE
    }
    if (is.vector(mods)) 
        mods <- cbind(mods)
    if (is.data.frame(mods)) 
        mods <- as.matrix(mods)
    if (is.character(mods)) 
        stop("Model matrix contains character variables.")
    if (!is.null(mods) && (nrow(mods) != k)) 
        stop("Number of rows of the model matrix does not match length of the outcome data.")
    if (very.verbose) 
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
        if (very.verbose) 
            message("Subsetting ...")
        mods <- mods[subset, , drop = FALSE]
        slab <- slab[subset]
        ids <- ids[subset]
    }
    if (anyDuplicated(slab)) 
        slab <- make.unique(as.character(slab))
    attr(yi, "slab") <- slab
    k <- length(yi)
    if (measure == "OR") {
        if (drop00) {
            id00 <- c(ai == 0L & ci == 0L) | c(bi == 0L & di == 
                0L)
            id00[is.na(id00)] <- FALSE
            ai[id00] <- NA
            bi[id00] <- NA
            ci[id00] <- NA
            di[id00] <- NA
        }
    }
    if (measure == "IRR") {
        if (drop00) {
            id00 <- c(x1i == 0L & x2i == 0L)
            id00[is.na(id00)] <- FALSE
            x1i[id00] <- NA
            x2i[id00] <- NA
        }
    }
    ai.f <- ai
    bi.f <- bi
    ci.f <- ci
    di.f <- di
    x1i.f <- x1i
    x2i.f <- x2i
    t1i.f <- t1i
    t2i.f <- t2i
    xi.f <- xi
    mi.f <- mi
    ti.f <- ti
    yi.f <- yi
    vi.f <- vi
    ni.f <- ni
    mods.f <- mods
    k.f <- k
    if (is.element(measure, c("OR"))) {
        aibicidimods.na <- is.na(ai) | is.na(bi) | is.na(ci) | 
            is.na(di) | if (is.null(mods)) 
            FALSE
        else apply(is.na(mods), 1, any)
        if (any(aibicidimods.na)) {
            if (very.verbose) 
                message("Handling NAs in table data ...")
            not.na <- !aibicidimods.na
            if (na.act == "na.omit" || na.act == "na.exclude" || 
                na.act == "na.pass") {
                ai <- ai[not.na]
                bi <- bi[not.na]
                ci <- ci[not.na]
                di <- di[not.na]
                mods <- mods[not.na, , drop = FALSE]
                k <- length(ai)
                warning("Studies with NAs omitted from model fitting.")
            }
            if (na.act == "na.fail") 
                stop("Missing values in studies.")
        }
        else {
            not.na <- rep(TRUE, k)
        }
    }
    if (is.element(measure, c("IRR"))) {
        x1ix2it1it2imods.na <- is.na(x1i) | is.na(x2i) | is.na(t1i) | 
            is.na(t2i) | if (is.null(mods)) 
            FALSE
        else apply(is.na(mods), 1, any)
        if (any(x1ix2it1it2imods.na)) {
            if (very.verbose) 
                message("Handling NAs in table data ...")
            not.na <- !x1ix2it1it2imods.na
            if (na.act == "na.omit" || na.act == "na.exclude" || 
                na.act == "na.pass") {
                x1i <- x1i[not.na]
                x2i <- x2i[not.na]
                t1i <- t1i[not.na]
                t2i <- t2i[not.na]
                mods <- mods[not.na, , drop = FALSE]
                k <- length(x1i)
                warning("Studies with NAs omitted from model fitting.")
            }
            if (na.act == "na.fail") 
                stop("Missing values in studies.")
        }
        else {
            not.na <- rep(TRUE, k)
        }
    }
    if (is.element(measure, c("PLO"))) {
        ximimods.na <- is.na(xi) | is.na(mi) | if (is.null(mods)) 
            FALSE
        else apply(is.na(mods), 1, any)
        if (any(ximimods.na)) {
            if (very.verbose) 
                message("Handling NAs in table data ...")
            not.na <- !ximimods.na
            if (na.act == "na.omit" || na.act == "na.exclude" || 
                na.act == "na.pass") {
                xi <- xi[not.na]
                mi <- mi[not.na]
                mods <- mods[not.na, , drop = FALSE]
                k <- length(xi)
                warning("Studies with NAs omitted from model fitting.")
            }
            if (na.act == "na.fail") 
                stop("Missing values in studies.")
        }
        else {
            not.na <- rep(TRUE, k)
        }
    }
    if (is.element(measure, c("IRLN"))) {
        xitimods.na <- is.na(xi) | is.na(ti) | if (is.null(mods)) 
            FALSE
        else apply(is.na(mods), 1, any)
        if (any(xitimods.na)) {
            if (very.verbose) 
                message("Handling NAs in table data ...")
            not.na <- !xitimods.na
            if (na.act == "na.omit" || na.act == "na.exclude" || 
                na.act == "na.pass") {
                xi <- xi[not.na]
                ti <- ti[not.na]
                mods <- mods[not.na, , drop = FALSE]
                k <- length(xi)
                warning("Studies with NAs omitted from model fitting.")
            }
            if (na.act == "na.fail") 
                stop("Missing values in studies.")
        }
        else {
            not.na <- rep(TRUE, k)
        }
    }
    if (k < 1) 
        stop("Processing terminated since k = 0.")
    mods.yi <- mods.f
    yivi.na <- is.na(yi) | is.na(vi) | if (is.null(mods.yi)) 
        FALSE
    else apply(is.na(mods.yi), 1, any)
    if (any(yivi.na)) {
        if (very.verbose) 
            message("Handling NAs in yi/vi ...")
        not.na.yivi <- !yivi.na
        if (na.act == "na.omit" || na.act == "na.exclude" || 
            na.act == "na.pass") {
            yi <- yi[not.na.yivi]
            ni <- ni[not.na.yivi]
            vi <- vi[not.na.yivi]
            mods.yi <- mods.f[not.na.yivi, , drop = FALSE]
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
    if (is.null(mods) && !intercept) {
        warning("Must either include an intercept and/or moderators in model.\n  Coerced intercept into the model.")
        intercept <- TRUE
    }
    if (intercept) {
        X <- cbind(intrcpt = rep(1, k), mods)
        X.f <- cbind(intrcpt = rep(1, k.f), mods.f)
        X.yi <- cbind(intrcpt = rep(1, k.yi), mods.yi)
    }
    else {
        X <- mods
        X.f <- mods.f
        X.yi <- mods.yi
    }
    is.int <- apply(X, 2, .is.int.func)
    if (any(is.int)) {
        int.incl <- TRUE
        int.indx <- which(is.int, arr.ind = TRUE)
        X <- cbind(intrcpt = 1, X[, -int.indx, drop = FALSE])
        X.f <- cbind(intrcpt = 1, X.f[, -int.indx, drop = FALSE])
        if (is.formula) 
            intercept <- TRUE
    }
    else {
        int.incl <- FALSE
    }
    is.int <- apply(X.yi, 2, .is.int.func)
    if (any(is.int)) {
        int.indx <- which(is.int, arr.ind = TRUE)
        X.yi <- cbind(intrcpt = 1, X.yi[, -int.indx, drop = FALSE])
    }
    tmp <- lm(rep(0, k) ~ X - 1)
    coef.na <- is.na(coef(tmp))
    if (any(coef.na)) {
        warning("Redundant predictors dropped from the model.")
        X <- X[, !coef.na, drop = FALSE]
        X.f <- X.f[, !coef.na, drop = FALSE]
    }
    tmp <- lm(yi ~ X.yi - 1)
    coef.na <- is.na(coef(tmp))
    if (any(coef.na)) 
        X.yi <- X.yi[, !coef.na, drop = FALSE]
    p <- NCOL(X)
    if ((p == 1L) && (all(sapply(X, identical, 1)))) {
        int.only <- TRUE
    }
    else {
        int.only <- FALSE
    }
    if (method == "FE" && p > k) 
        stop("Number of parameters to be estimated is larger than the number of observations.")
    if (method != "FE" && (p + 1) > k) 
        stop("Number of parameters to be estimated is larger than the number of observations.")
    if (is.null(btt)) {
        if (p > 1) {
            if (int.incl) {
                btt <- seq.int(from = 2, to = p)
            }
            else {
                btt <- seq_len(p)
            }
        }
        else {
            btt <- 1
        }
    }
    else {
        btt <- btt[(btt >= 1) & (btt <= p)]
        btt <- unique(round(btt))
        if (length(btt) == 0L) 
            stop("Non-existent coefficients specified via 'btt'.")
    }
    m <- length(btt)
    con <- list(verbose = FALSE, optimizer = "optim", optmethod = "BFGS", 
        scale = TRUE, tol = 1e-07, dnchgcalc = "dFNCHypergeo", 
        dnchgprec = 1e-10)
    con.pos <- pmatch(names(control), names(con))
    con[c(na.omit(con.pos))] <- control[!is.na(con.pos)]
    if (verbose) 
        con$verbose <- verbose
    verbose <- con$verbose
    pos.optCtrl <- pmatch(names(control), "optCtrl", nomatch = 0)
    if (sum(pos.optCtrl) > 0) {
        optCtrl <- control[[which(pos.optCtrl == 1)]]
    }
    else {
        optCtrl <- list()
    }
    if (con$optimizer == "optim") {
        con.pos <- pmatch(names(optCtrl), "REPORT", nomatch = 0)
        if (sum(con.pos) > 0) {
            optCtrl[which(con.pos == 1)] <- 1
            names(optCtrl)[which(con.pos == 1)] <- "REPORT"
        }
        else {
            optCtrl$REPORT <- 1
        }
        optCtrl$trace <- con$verbose
    }
    if (con$optimizer == "nlminb") 
        optCtrl$trace <- ifelse(con$verbose > 0, 1, 0)
    if (is.element(con$optimizer, c("uobyqa", "newuoa", "bobyqa"))) 
        optCtrl$iprint <- ifelse(con$verbose > 0, 3, 0)
    pos.clogitCtrl <- pmatch(names(control), "clogitCtrl", nomatch = 0)
    if (sum(pos.clogitCtrl) > 0) {
        clogitCtrl <- control[[which(pos.clogitCtrl == 1)]]
    }
    else {
        clogitCtrl <- list()
    }
    pos.clogisticCtrl <- pmatch(names(control), "clogisticCtrl", 
        nomatch = 0)
    if (sum(pos.clogisticCtrl) > 0) {
        clogisticCtrl <- control[[which(pos.clogisticCtrl == 
            1)]]
    }
    else {
        clogisticCtrl <- list()
    }
    pos.glmCtrl <- pmatch(names(control), "glmCtrl", nomatch = 0)
    if (sum(pos.glmCtrl) > 0) {
        glmCtrl <- control[[which(pos.glmCtrl == 1)]]
    }
    else {
        glmCtrl <- list()
    }
    glmCtrl$trace <- ifelse(con$verbose > 0, TRUE, FALSE)
    pos.glmerCtrl <- pmatch(names(control), "glmerCtrl", nomatch = 0)
    if (sum(pos.glmerCtrl) > 0) {
        glmerCtrl <- control[[which(pos.glmerCtrl == 1)]]
    }
    else {
        glmerCtrl <- list()
    }
    pos.intCtrl <- pmatch(names(control), "intCtrl", nomatch = 0)
    if (sum(pos.intCtrl) > 0) {
        intCtrl <- control[[which(pos.intCtrl == 1)]]
    }
    else {
        intCtrl <- list()
    }
    con.pos <- pmatch(names(intCtrl), "lower", nomatch = 0)
    if (sum(con.pos) > 0) {
        intCtrl[which(con.pos == 1)] <- -Inf
        names(intCtrl)[which(con.pos == 1)] <- "lower"
    }
    else {
        intCtrl$lower <- -Inf
    }
    con.pos <- pmatch(names(intCtrl), "upper", nomatch = 0)
    if (sum(con.pos) > 0) {
        intCtrl[which(con.pos == 1)] <- Inf
        names(intCtrl)[which(con.pos == 1)] <- "upper"
    }
    else {
        intCtrl$upper <- Inf
    }
    con.pos <- pmatch(names(intCtrl), "subdivisions", nomatch = 0)
    if (sum(con.pos) > 0) {
        intCtrl[which(con.pos == 1)] <- 100L
        names(intCtrl)[which(con.pos == 1)] <- "subdivisions"
    }
    else {
        intCtrl$subdivisions <- 100L
    }
    con.pos <- pmatch(names(intCtrl), "rel.tol", nomatch = 0)
    if (sum(con.pos) > 0) {
        intCtrl[which(con.pos == 1)] <- .Machine$double.eps^0.25
        names(intCtrl)[which(con.pos == 1)] <- "rel.tol"
    }
    else {
        intCtrl$rel.tol <- .Machine$double.eps^0.25
    }
    if (!is.element(con$optimizer, c("optim", "nlminb", "uobyqa", 
        "newuoa", "bobyqa", "clogit", "clogistic"))) 
        stop("Unknown optimizer specified.")
    if (con$dnchgcalc != "dnoncenhypergeom" && con$dnchgcalc != 
        "dFNCHypergeo") 
        stop("Unknown dnchgcalc method specified.")
    if (is.element(con$optimizer, c("clogit", "clogistic")) && 
        method == "ML") 
        stop("Cannot use 'clogit' or 'clogistic' with method='ML'.")
    if (is.element(measure, c("OR", "IRR"))) {
        if ((model == "UM.FS" && method == "ML") || (model == 
            "UM.RS") || (model == "CM.AL" && method == "ML") || 
            (model == "CM.EL" && method == "ML")) {
            if (!requireNamespace("lme4", quietly = TRUE)) 
                stop("Please install the 'lme4' package to fit this model.")
        }
    }
    if (is.element(measure, c("PLO", "IRLN")) && method == "ML") {
        if (!requireNamespace("lme4", quietly = TRUE)) 
            stop("Please install the 'lme4' package to fit this model.")
    }
    if (measure == "OR" && model == "CM.EL") {
        if (is.element(con$optimizer, c("uobyqa", "newuoa", "bobyqa"))) {
            if (!requireNamespace("minqa", quietly = TRUE)) 
                stop("Please install the 'minqa' package to fit this model.")
            minqa <- get(con$optimizer, envir = loadNamespace("minqa"))
            con$optimizer <- "minqa"
        }
        if (con$optimizer == "optim" || con$optimizer == "nlminb" || 
            con$optimizer == "minqa") {
            if (!requireNamespace("numDeriv", quietly = TRUE)) 
                stop("Please install the 'numDeriv' package to fit this model.")
            if (con$dnchgcalc == "dFNCHypergeo") {
                if (!requireNamespace("BiasedUrn", quietly = TRUE)) 
                  stop("Please install the 'BiasedUrn' package to fit this model.")
            }
        }
        if (con$optimizer == "clogit") {
            if (!requireNamespace("survival", quietly = TRUE)) 
                stop("Please install the 'survival' package to fit this model.")
            coxph <- survival::coxph
            Surv <- survival::Surv
        }
        if (con$optimizer == "clogistic") {
            if (!requireNamespace("Epi", quietly = TRUE)) 
                stop("Please install the 'Epi' package to fit this model.")
        }
    }
    if (any(eigen(crossprod(X), symmetric = TRUE, only.values = TRUE)$values <= 
        con$tol)) 
        stop("Model matrix not of full rank. Cannot fit model.")
    se.tau2 <- I2 <- H2 <- QE <- QEp <- NA
    alpha <- ifelse(level > 1, (100 - level)/100, 1 - level)
    if (!int.only && int.incl && con$scale) {
        Xsave <- X
        meanX <- colMeans(X[, 2:p, drop = FALSE])
        sdX <- apply(X[, 2:p, drop = FALSE], 2, sd)
        is.d <- apply(X, 2, function(x) all(sapply(x, identical, 
            0) | sapply(x, identical, 1)))
        X[, !is.d] <- apply(X[, !is.d, drop = FALSE], 2, scale)
    }
    if (is.element(measure, c("OR", "IRR"))) {
        if (is.element(model, c("UM.FS", "UM.RS"))) {
            if (measure == "OR") {
                dat.grp <- cbind(xi = c(rbind(ai, ci)), mi = c(rbind(bi, 
                  di)))
                dat.fam <- binomial
                dat.off <- NULL
            }
            if (measure == "IRR") {
                dat.grp <- cbind(xi = c(rbind(x1i, x2i)))
                dat.fam <- poisson
                dat.off <- log(c(rbind(t1i, t2i)))
            }
            group1 <- rep(c(1, 0), times = k)
            group2 <- rep(c(0, 1), times = k)
            group12 <- rep(c(1/2, -1/2), times = k)
            study <- factor(rep(seq_len(k), each = 2))
            const <- cbind(rep(1, 2 * k))
            X.fit <- X[rep(seq(k), each = 2), , drop = FALSE]
            X.fit <- cbind(group1 * X.fit[, , drop = FALSE])
            row.names(X.fit) <- seq_len(2 * k)
            if (model == "UM.FS") {
                if (verbose) 
                  message("Fitting FE model ...")
                if (k > 1) {
                  res.FE <- try(glm(dat.grp ~ -1 + X.fit + study, 
                    offset = dat.off, family = dat.fam, control = glmCtrl), 
                    silent = !verbose)
                }
                else {
                  res.FE <- try(glm(dat.grp ~ -1 + X.fit + const, 
                    offset = dat.off, family = dat.fam, control = glmCtrl), 
                    silent = !verbose)
                }
                if (inherits(res.FE, "try-error")) 
                  stop("Cannot fit FE model.")
                if (verbose) 
                  message("Fitting saturated model ...")
                if (k > 1) {
                  X.QE <- model.matrix(~-1 + X.fit + study + 
                    study:group1)
                  res.QE <- try(glm(dat.grp ~ -1 + X.QE, offset = dat.off, 
                    family = dat.fam, control = glmCtrl), silent = !verbose)
                }
                else {
                  res.QE <- res.FE
                }
                if (inherits(res.QE, "try-error")) 
                  stop("Cannot fit saturated model.")
                ll.FE <- c(logLik(res.FE))
                ll.QE <- c(logLik(res.QE))
                b2.QE <- cbind(na.omit(coef(res.QE)[-seq_len(k + 
                  p)]))
                vb2.QE <- vcov(res.QE)[-seq_len(k + p), -seq_len(k + 
                  p), drop = FALSE]
                if (method == "ML") {
                  if (verbose) 
                    message("Fitting ML model ...")
                  res.ML <- try(lme4::glmer(dat.grp ~ -1 + X.fit + 
                    study + (group12 - 1 | study), offset = dat.off, 
                    family = dat.fam, nAGQ = nAGQ, verbose = verbose, 
                    control = do.call(lme4::glmerControl, glmerCtrl)), 
                    silent = !verbose)
                  if (inherits(res.ML, "try-error")) 
                    stop("Cannot fit ML model.")
                  ll.ML <- ll.QE - 1/2 * deviance(res.ML)
                }
                if (method == "FE") {
                  b <- cbind(coef(res.FE)[seq_len(p)])
                  vb <- vcov(res.FE)[seq_len(p), seq_len(p), 
                    drop = FALSE]
                  tau2 <- 0
                  sigma2 <- NA
                  parms <- p + k
                  p.eff <- p + k
                  k.eff <- 2 * k
                }
                if (method == "ML") {
                  b <- cbind(lme4::fixef(res.ML)[seq_len(p)])
                  vb <- as.matrix(vcov(res.ML))[seq_len(p), seq_len(p), 
                    drop = FALSE]
                  tau2 <- lme4::VarCorr(res.ML)[[1]][1]
                  sigma2 <- NA
                  parms <- p + k + 1
                  p.eff <- p + k
                  k.eff <- 2 * k
                }
            }
            if (model == "UM.RS") {
                if (verbose) 
                  message("Fitting FE model ...")
                res.FE <- try(lme4::glmer(dat.grp ~ -1 + X.fit + 
                  const + (1 | study), offset = dat.off, family = dat.fam, 
                  nAGQ = nAGQ, verbose = verbose, control = do.call(lme4::glmerControl, 
                    glmerCtrl)), silent = !verbose)
                if (inherits(res.FE, "try-error")) 
                  stop("Cannot fit FE model.")
                if (verbose) 
                  message("Fitting saturated model ...")
                if (k > 1) {
                  X.QE <- model.matrix(~-1 + X.fit + const + 
                    study:group1)
                  res.QE <- try(glm(dat.grp ~ -1 + X.QE, offset = dat.off, 
                    family = dat.fam, control = glmCtrl), silent = !verbose)
                  X.QE <- X.QE[, !is.na(coef(res.QE)), drop = FALSE]
                  res.QE <- try(lme4::glmer(dat.grp ~ -1 + X.QE + 
                    (1 | study), offset = dat.off, family = dat.fam, 
                    start = c(sqrt(lme4::VarCorr(res.FE)[[1]][1])), 
                    nAGQ = nAGQ, verbose = verbose, control = do.call(lme4::glmerControl, 
                      glmerCtrl)), silent = !verbose)
                }
                else {
                  res.QE <- res.FE
                }
                if (inherits(res.QE, "try-error")) 
                  stop("Cannot fit saturated model.")
                ll.FE <- c(logLik(res.FE))
                ll.QE <- c(logLik(res.QE))
                b2.QE <- cbind(lme4::fixef(res.QE)[-seq_len(p + 
                  1)])
                vb2.QE <- as.matrix(vcov(res.QE))[-seq_len(p + 
                  1), -seq_len(p + 1), drop = FALSE]
                if (method == "ML") {
                  if (verbose) 
                    message("Fitting ML model ...")
                  res.ML <- try(lme4::glmer(dat.grp ~ -1 + X.fit + 
                    const + (1 | study) + (group12 - 1 | study), 
                    offset = dat.off, family = dat.fam, nAGQ = nAGQ, 
                    verbose = verbose, control = do.call(lme4::glmerControl, 
                      glmerCtrl)), silent = !verbose)
                  if (inherits(res.ML, "try-error")) 
                    stop("Cannot fit ML model.")
                  ll.ML <- c(logLik(res.ML))
                }
                if (method == "FE") {
                  b <- cbind(lme4::fixef(res.FE)[seq_len(p)])
                  vb <- as.matrix(vcov(res.FE))[seq_len(p), seq_len(p), 
                    drop = FALSE]
                  tau2 <- 0
                  sigma2 <- lme4::VarCorr(res.FE)[[1]][1]
                  parms <- p + 1 + 1
                  p.eff <- p + 1
                  k.eff <- 2 * k
                }
                if (method == "ML") {
                  b <- cbind(lme4::fixef(res.ML)[seq_len(p)])
                  vb <- as.matrix(vcov(res.ML))[seq_len(p), seq_len(p), 
                    drop = FALSE]
                  tau2 <- lme4::VarCorr(res.ML)[[2]][1]
                  sigma2 <- lme4::VarCorr(res.ML)[[1]][1]
                  parms <- p + 1 + 2
                  p.eff <- p + 1
                  k.eff <- 2 * k
                }
            }
        }
        if ((measure == "IRR" && model == "CM.EL") || (measure == 
            "OR" && model == "CM.AL") || (measure == "OR" && 
            model == "CM.EL")) {
            if (measure == "OR") {
                dat.grp <- cbind(xi = ai, mi = ci)
                dat.off <- log((ai + bi)/(ci + di))
            }
            if (measure == "IRR") {
                dat.grp <- cbind(xi = x1i, mi = x2i)
                dat.off <- log(t1i/t2i)
            }
            study <- factor(seq_len(k))
            X.fit <- X
            if (verbose) 
                message("Fitting FE model ...")
            res.FE <- try(glm(dat.grp ~ -1 + X.fit, offset = dat.off, 
                family = binomial, control = glmCtrl), silent = !verbose)
            if (inherits(res.FE, "try-error")) 
                stop("Cannot fit FE model.")
            if (verbose) 
                message("Fitting saturated model ...")
            if (k > 1) {
                X.QE <- model.matrix(~-1 + X.fit + study)
                res.QE <- try(glm(dat.grp ~ -1 + X.QE, offset = dat.off, 
                  family = binomial, control = glmCtrl), silent = !verbose)
            }
            else {
                res.QE <- res.FE
            }
            if (inherits(res.QE, "try-error")) 
                stop("Cannot fit saturated model.")
            ll.FE <- c(logLik(res.FE))
            ll.QE <- c(logLik(res.QE))
            b2.QE <- cbind(na.omit(coef(res.QE)[-seq_len(p)]))
            vb2.QE <- vcov(res.QE)[-seq_len(p), -seq_len(p), 
                drop = FALSE]
            if (method == "ML") {
                if (verbose) 
                  message("Fitting ML model ...")
                if (verbose) {
                  res.ML <- try(lme4::glmer(dat.grp ~ -1 + X.fit + 
                    (1 | study), offset = dat.off, family = binomial, 
                    nAGQ = nAGQ, verbose = verbose, control = do.call(lme4::glmerControl, 
                      glmerCtrl)), silent = !verbose)
                }
                else {
                  res.ML <- suppressMessages(try(lme4::glmer(dat.grp ~ 
                    -1 + X.fit + (1 | study), offset = dat.off, 
                    family = binomial, nAGQ = nAGQ, verbose = verbose, 
                    control = do.call(lme4::glmerControl, glmerCtrl)), 
                    silent = !verbose))
                }
                if (inherits(res.ML, "try-error")) 
                  stop("Cannot fit ML model.")
                ll.ML <- ll.QE - 1/2 * deviance(res.ML)
            }
            if (method == "FE") {
                b <- cbind(coef(res.FE)[seq_len(p)])
                vb <- vcov(res.FE)[seq_len(p), seq_len(p), drop = FALSE]
                tau2 <- 0
                sigma2 <- NA
                parms <- p
                p.eff <- p
                k.eff <- k
            }
            if (method == "ML") {
                b <- cbind(lme4::fixef(res.ML)[seq_len(p)])
                vb <- as.matrix(vcov(res.ML))[seq_len(p), seq_len(p), 
                  drop = FALSE]
                tau2 <- lme4::VarCorr(res.ML)[[1]][1]
                sigma2 <- NA
                parms <- p + 1
                p.eff <- p
                k.eff <- k
            }
        }
        if (measure == "OR" && model == "CM.EL") {
            if (verbose) 
                message("Fitting FE model ...")
            if (con$optimizer == "optim" || con$optimizer == 
                "nlminb" || con$optimizer == "minqa") {
                if (con$optimizer == "optim") {
                  res.FE <- try(optim(par = c(coef(res.FE)[seq_len(p)], 
                    0), .dnchg, method = con$optmethod, hessian = TRUE, 
                    ai = ai, bi = bi, ci = ci, di = di, X.fit = X.fit, 
                    random = FALSE, verbose = verbose, digits = digits, 
                    dnchgcalc = con$dnchgcalc, dnchgprec = con$dnchgprec, 
                    control = optCtrl), silent = !verbose)
                }
                if (con$optimizer == "nlminb") {
                  res.FE <- try(nlminb(start = c(coef(res.FE)[seq_len(p)], 
                    0), .dnchg, ai = ai, bi = bi, ci = ci, di = di, 
                    X.fit = X.fit, random = FALSE, verbose = verbose, 
                    digits = digits, dnchgcalc = con$dnchgcalc, 
                    dnchgprec = con$dnchgprec, control = optCtrl), 
                    silent = !verbose)
                }
                if (con$optimizer == "minqa") {
                  res.FE <- try(minqa(par = c(coef(res.FE)[seq_len(p)], 
                    0), .dnchg, ai = ai, bi = bi, ci = ci, di = di, 
                    X.fit = X.fit, random = FALSE, verbose = verbose, 
                    digits = digits, dnchgcalc = con$dnchgcalc, 
                    dnchgprec = con$dnchgprec, control = optCtrl), 
                    silent = !verbose)
                }
                if (con$optimizer == "optim" || con$optimizer == 
                  "nlminb") {
                  if (inherits(res.FE, "try-error") || res.FE$convergence != 
                    0) 
                    stop("Cannot fit FE model.")
                }
                if (con$optimizer == "minqa") {
                  if (inherits(res.FE, "try-error") || res.FE$ierr != 
                    0) 
                    stop("Cannot fit FE model.")
                }
                if (very.verbose) 
                  message("Computing Hessian ...")
                h.FE <- numDeriv::hessian(.dnchg, x = res.FE$par, 
                  ai = ai, bi = bi, ci = ci, di = di, X.fit = X.fit, 
                  random = FALSE, verbose = verbose, digits = digits, 
                  dnchgcalc = con$dnchgcalc, dnchgprec = con$dnchgprec)
                if (verbose) 
                  message("Fitting saturated model ...")
                QEconv <- TRUE
                QE.Wld <- NA
                if (k > 1) {
                  is.aliased <- is.na(coef(res.QE))
                  X.QE <- X.QE[, !is.aliased, drop = FALSE]
                  if (con$optimizer == "optim") {
                    res.QE <- try(optim(par = c(coef(res.QE)[!is.aliased], 
                      0), .dnchg, method = con$optmethod, hessian = TRUE, 
                      ai = ai, bi = bi, ci = ci, di = di, X.fit = X.QE, 
                      random = FALSE, verbose = verbose, digits = digits, 
                      dnchgcalc = con$dnchgcalc, dnchgprec = con$dnchgprec, 
                      control = optCtrl), silent = !verbose)
                  }
                  if (con$optimizer == "nlminb") {
                    res.QE <- try(nlminb(start = c(coef(res.QE)[!is.aliased], 
                      0), .dnchg, ai = ai, bi = bi, ci = ci, 
                      di = di, X.fit = X.QE, random = FALSE, 
                      verbose = verbose, digits = digits, dnchgcalc = con$dnchgcalc, 
                      dnchgprec = con$dnchgprec, control = optCtrl), 
                      silent = !verbose)
                  }
                  if (con$optimizer == "minqa") {
                    res.QE <- try(minqa(par = c(coef(res.QE)[!is.aliased], 
                      0), .dnchg, ai = ai, bi = bi, ci = ci, 
                      di = di, X.fit = X.QE, random = FALSE, 
                      verbose = verbose, digits = digits, dnchgcalc = con$dnchgcalc, 
                      dnchgprec = con$dnchgprec, control = optCtrl), 
                      silent = !verbose)
                  }
                  if (con$optimizer == "optim" || con$optimizer == 
                    "nlminb") {
                    if (inherits(res.QE, "try-error") || res.QE$convergence != 
                      0) {
                      warning("Cannot fit saturated model.")
                      QEconv <- FALSE
                      res.QE <- NULL
                    }
                  }
                  if (con$optimizer == "minqa") {
                    if (inherits(res.QE, "try-error") || res.QE$ierr != 
                      0) {
                      warning("Cannot fit saturated model.")
                      QEconv <- FALSE
                      res.QE <- NULL
                    }
                  }
                  if (QEconv) {
                    if (very.verbose) 
                      message("Computing Hessian ...")
                    h.QE <- numDeriv::hessian(.dnchg, x = res.QE$par, 
                      ai = ai, bi = bi, ci = ci, di = di, X.fit = X.QE, 
                      random = FALSE, verbose = verbose, digits = digits, 
                      dnchgcalc = con$dnchgcalc, dnchgprec = con$dnchgprec)
                  }
                }
                else {
                  res.QE <- res.FE
                  h.QE <- h.FE
                }
                if (con$optimizer == "optim") {
                  ll.FE <- -1 * res.FE$value
                  ll.QE <- -1 * res.QE$value
                }
                if (con$optimizer == "nlminb") {
                  ll.FE <- -1 * res.FE$objective
                  ll.QE <- -1 * res.QE$objective
                }
                if (con$optimizer == "minqa") {
                  ll.FE <- -1 * res.FE$fval
                  ll.QE <- -1 * res.QE$fval
                }
                if (!QEconv) 
                  ll.QE <- NA
                if (QEconv) {
                  b2.QE <- res.QE$par
                  hessian <- h.QE
                  p.QE <- length(b2.QE)
                  b2.QE <- b2.QE[-p.QE]
                  hessian <- hessian[-p.QE, -p.QE, drop = FALSE]
                  p.QE <- length(b2.QE)
                  is.0 <- colSums(hessian == 0L) == p.QE
                  b2.QE <- b2.QE[!is.0]
                  hessian <- hessian[!is.0, !is.0, drop = FALSE]
                  b2.QE <- cbind(b2.QE[-seq_len(p)])
                  h.A <- hessian[seq_len(p), seq_len(p), drop = FALSE]
                  h.B <- hessian[seq_len(p), -seq_len(p), drop = FALSE]
                  h.C <- hessian[-seq_len(p), seq_len(p), drop = FALSE]
                  h.D <- hessian[-seq_len(p), -seq_len(p), drop = FALSE]
                  chol.h.A <- try(chol(h.A), silent = !verbose)
                  if (class(chol.h.A) == "try-error") {
                    warning("Cannot invert Hessian for saturated model.")
                  }
                  else {
                    Ivb2.QE <- h.D - h.C %*% chol2inv(chol.h.A) %*% 
                      h.B
                    QE.Wld <- c(t(b2.QE) %*% Ivb2.QE %*% b2.QE)
                  }
                }
            }
            if (con$optimizer == "clogit" || con$optimizer == 
                "clogistic") {
                event <- unlist(lapply(seq_len(k), function(i) c(rep.int(1, 
                  ai[i]), rep.int(0, bi[i]), rep.int(1, ci[i]), 
                  rep.int(0, di[i]))))
                group1 <- unlist(lapply(seq_len(k), function(i) c(rep.int(1, 
                  ai[i]), rep.int(1, bi[i]), rep.int(0, ci[i]), 
                  rep.int(0, di[i]))))
                study.l <- factor(rep(seq_len(k), times = ni))
                X.fit.l <- X[rep(seq_len(k), times = ni), , drop = FALSE]
                X.fit.l <- cbind(group1 * X.fit.l)
                const <- rep(1, length(event))
                if (k > 1) {
                  if (con$optimizer == "clogit") {
                    args.clogit <- clogitCtrl
                    args.clogit$formula <- event ~ X.fit.l + 
                      strata(study.l)
                    res.FE <- try(do.call(survival::clogit, args.clogit), 
                      silent = !verbose)
                  }
                  if (con$optimizer == "clogistic") {
                    args.clogistic <- clogisticCtrl
                    args.clogistic$formula <- event ~ X.fit.l
                    args.clogistic$strata <- study.l
                    res.FE <- try(do.call(Epi::clogistic, args.clogistic), 
                      silent = !verbose)
                  }
                }
                else {
                  stop(paste("Cannot use '", con$optimizer, "' optimizer when k=1.", 
                    sep = ""))
                }
                if (inherits(res.FE, "try-error")) 
                  stop("Cannot fit FE model.")
                if (verbose) 
                  message("Fitting saturated model ...")
                X.QE.l <- model.matrix(~-1 + X.fit.l + study.l:group1)
                X.QE.l <- X.QE.l[, !is.na(coef(res.QE)), drop = FALSE]
                X.QE <- X.QE[, !is.na(coef(res.QE)), drop = FALSE]
                if (con$optimizer == "clogit") {
                  args.clogit <- clogitCtrl
                  args.clogit$formula <- event ~ X.QE.l + strata(study.l)
                  if (verbose) {
                    res.QE <- try(do.call(survival::clogit, args.clogit), 
                      silent = !verbose)
                  }
                  else {
                    res.QE <- suppressWarnings(try(do.call(survival::clogit, 
                      args.clogit), silent = !verbose))
                  }
                }
                if (con$optimizer == "clogistic") {
                  args.clogistic <- clogisticCtrl
                  args.clogistic$formula <- event ~ X.QE.l
                  args.clogistic$strata <- study.l
                  res.QE <- try(do.call(Epi::clogistic, args.clogistic), 
                    silent = !verbose)
                }
                if (inherits(res.QE, "try-error")) 
                  stop("Cannot fit saturated model.")
                ll.FE <- -1 * .dnchg(c(cbind(coef(res.FE)), 0), 
                  ai = ai, bi = bi, ci = ci, di = di, X.fit = X.fit, 
                  random = FALSE, verbose = verbose, digits = digits, 
                  dnchgcalc = con$dnchgcalc, dnchgprec = con$dnchgprec)
                ll.QE <- -1 * .dnchg(c(cbind(coef(res.QE)), 0), 
                  ai = ai, bi = bi, ci = ci, di = di, X.fit = X.QE, 
                  random = FALSE, verbose = verbose, digits = digits, 
                  dnchgcalc = con$dnchgcalc, dnchgprec = con$dnchgprec)
                b2.QE <- cbind(coef(res.QE)[-seq_len(p)])
                vb2.QE <- vcov(res.QE)[-seq_len(p), -seq_len(p), 
                  drop = FALSE]
            }
            if (method == "ML") {
                if (verbose) 
                  message("Fitting ML model ...")
                if (con$optimizer == "optim") {
                  res.ML <- try(optim(par = c(b, log(tau2 + 0.001)), 
                    .dnchg, method = con$optmethod, hessian = TRUE, 
                    ai = ai, bi = bi, ci = ci, di = di, X.fit = X.fit, 
                    random = TRUE, verbose = verbose, digits = digits, 
                    dnchgcalc = con$dnchgcalc, dnchgprec = con$dnchgprec, 
                    intCtrl = intCtrl, control = optCtrl), silent = !verbose)
                }
                if (con$optimizer == "nlminb") {
                  res.ML <- try(nlminb(start = c(b, log(tau2 + 
                    0.001)), .dnchg, ai = ai, bi = bi, ci = ci, 
                    di = di, X.fit = X.fit, random = TRUE, verbose = verbose, 
                    digits = digits, dnchgcalc = con$dnchgcalc, 
                    dnchgprec = con$dnchgprec, intCtrl = intCtrl, 
                    control = optCtrl), silent = !verbose)
                }
                if (con$optimizer == "minqa") {
                  res.ML <- try(minqa(par = c(b, log(tau2 + 0.001)), 
                    .dnchg, ai = ai, bi = bi, ci = ci, di = di, 
                    X.fit = X.fit, random = TRUE, verbose = verbose, 
                    digits = digits, dnchgcalc = con$dnchgcalc, 
                    dnchgprec = con$dnchgprec, intCtrl = intCtrl, 
                    control = optCtrl), silent = !verbose)
                }
                if (con$optimizer == "optim" || con$optimizer == 
                  "nlminb") {
                  if (inherits(res.ML, "try-error") || res.ML$convergence != 
                    0) 
                    stop("Cannot fit ML model.")
                }
                if (con$optimizer == "minqa") {
                  if (inherits(res.ML, "try-error") || res.ML$ierr != 
                    0) 
                    stop("Cannot fit ML model.")
                }
                if (very.verbose) 
                  message("Computing Hessian ...")
                h.ML <- numDeriv::hessian(.dnchg, x = res.ML$par, 
                  method.args = list(r = 8), ai = ai, bi = bi, 
                  ci = ci, di = di, X.fit = X.fit, random = TRUE, 
                  verbose = verbose, digits = digits, dnchgcalc = con$dnchgcalc, 
                  dnchgprec = con$dnchgprec, intCtrl = intCtrl)
                if (con$optimizer == "optim") {
                  ll.ML <- -1 * res.ML$value
                }
                if (con$optimizer == "nlminb") {
                  ll.ML <- -1 * res.ML$objective
                }
                if (con$optimizer == "minqa") {
                  ll.ML <- -1 * res.ML$fval
                }
            }
            if (method == "FE") {
                if (con$optimizer == "optim" || con$optimizer == 
                  "nlminb" || con$optimizer == "minqa") {
                  b <- cbind(res.FE$par[seq_len(p)])
                  chol.h <- try(chol(h.FE[seq_len(p), seq_len(p)]), 
                    silent = !verbose)
                  if (class(chol.h) == "try-error") {
                    warning("Choleski factorization of Hessian failed. Trying inversion via QR decomposition.")
                    vb <- try(qr.solve(h.FE[seq_len(p), seq_len(p)]), 
                      silent = !verbose)
                    if (class(vb) == "try-error") 
                      stop("Cannot invert Hessian for ML model.")
                  }
                  else {
                    vb <- chol2inv(chol.h)
                  }
                }
                if (con$optimizer == "clogit" || con$optimizer == 
                  "clogistic") {
                  b <- cbind(coef(res.FE)[seq_len(p)])
                  vb <- vcov(res.FE)[seq_len(p), seq_len(p), 
                    drop = FALSE]
                }
                tau2 <- 0
                sigma2 <- NA
                parms <- p
                p.eff <- p
                k.eff <- k
            }
            if (method == "ML") {
                b <- cbind(res.ML$par[seq_len(p)])
                chol.h <- try(chol(h.ML), silent = !verbose)
                if (class(chol.h) == "try-error") {
                  warning("Choleski factorization of Hessian failed. Trying inversion via QR decomposition.")
                  vb.f <- try(qr.solve(h.ML), silent = !verbose)
                  if (class(vb.f) == "try-error") 
                    stop("Cannot invert Hessian for ML model.")
                }
                else {
                  vb.f <- chol2inv(chol.h)
                }
                vb <- vb.f[seq_len(p), seq_len(p), drop = FALSE]
                tau2 <- exp(res.ML$par[p + 1])
                sigma2 <- NA
                parms <- p + 1
                p.eff <- p
                k.eff <- k
                if (vb.f[p + 1, p + 1] >= 0) {
                  se.tau2 <- sqrt(vb.f[p + 1, p + 1]) * tau2
                }
                else {
                  se.tau2 <- NA
                }
            }
        }
    }
    if (is.element(measure, c("PLO", "IRLN"))) {
        if (measure == "PLO") {
            dat.grp <- cbind(xi = xi, mi = mi)
            dat.fam <- binomial
            dat.off <- NULL
        }
        if (measure == "IRLN") {
            dat.grp <- xi
            dat.fam <- poisson
            dat.off <- log(ti)
        }
        study <- factor(seq_len(k))
        X.fit <- X
        if (verbose) 
            message("Fitting FE model ...")
        res.FE <- try(glm(dat.grp ~ -1 + X.fit, offset = dat.off, 
            family = dat.fam, control = glmCtrl), silent = !verbose)
        if (inherits(res.FE, "try-error")) 
            stop("Cannot fit FE model.")
        if (verbose) 
            message("Fitting saturated model ...")
        if (k > 1) {
            X.QE <- model.matrix(~-1 + X.fit + study)
            if (verbose) {
                res.QE <- try(glm(dat.grp ~ -1 + X.QE, offset = dat.off, 
                  family = dat.fam, control = glmCtrl), silent = !verbose)
            }
            else {
                res.QE <- suppressWarnings(try(glm(dat.grp ~ 
                  -1 + X.QE, offset = dat.off, family = dat.fam, 
                  control = glmCtrl), silent = !verbose))
            }
        }
        else {
            res.QE <- res.FE
        }
        if (inherits(res.QE, "try-error")) 
            stop("Cannot fit saturated model.")
        ll.FE <- c(logLik(res.FE))
        ll.QE <- c(logLik(res.QE))
        b2.QE <- cbind(na.omit(coef(res.QE)[-seq_len(p)]))
        vb2.QE <- vcov(res.QE)[-seq_len(p), -seq_len(p), drop = FALSE]
        if (method == "ML") {
            if (verbose) 
                message("Fitting ML model ...")
            if (verbose) {
                res.ML <- try(lme4::glmer(dat.grp ~ -1 + X.fit + 
                  (1 | study), offset = dat.off, family = dat.fam, 
                  nAGQ = nAGQ, verbose = verbose, control = do.call(lme4::glmerControl, 
                    glmerCtrl)), silent = !verbose)
            }
            else {
                res.ML <- suppressMessages(try(lme4::glmer(dat.grp ~ 
                  -1 + X.fit + (1 | study), offset = dat.off, 
                  family = dat.fam, nAGQ = nAGQ, verbose = verbose, 
                  control = do.call(lme4::glmerControl, glmerCtrl)), 
                  silent = !verbose))
            }
            if (inherits(res.ML, "try-error")) 
                stop("Cannot fit ML model.")
            ll.ML <- ll.QE - 1/2 * deviance(res.ML)
        }
        if (method == "FE") {
            b <- cbind(coef(res.FE)[seq_len(p)])
            vb <- vcov(res.FE)[seq_len(p), seq_len(p), drop = FALSE]
            tau2 <- 0
            sigma2 <- NA
            parms <- p
            p.eff <- p
            k.eff <- k
        }
        if (method == "ML") {
            b <- cbind(lme4::fixef(res.ML)[seq_len(p)])
            vb <- as.matrix(vcov(res.ML))[seq_len(p), seq_len(p), 
                drop = FALSE]
            tau2 <- lme4::VarCorr(res.ML)[[1]][1]
            sigma2 <- NA
            parms <- p + 1
            p.eff <- p
            k.eff <- k
        }
    }
    if (very.verbose) 
        message("Heterogeneity testing ...")
    if (measure != "OR" || model != "CM.EL" || con$optimizer != 
        "optim") {
        if (dim(vb2.QE)[1] > 0) {
            chol.h <- try(chol(vb2.QE), silent = !verbose)
            if (class(chol.h) == "try-error") {
                warning("Cannot invert Hessian for saturated model.")
                QE.Wld <- NA
            }
            else {
                QE.Wld <- c(t(b2.QE) %*% chol2inv(chol.h) %*% 
                  b2.QE)
            }
        }
        else {
            QE.Wld <- 0
        }
    }
    QE.LRT <- -2 * (ll.FE - ll.QE)
    QE.Wld[QE.Wld < 0] <- 0
    QE.LRT[QE.LRT < 0] <- 0
    QE.df <- k - p
    if (QE.df > 0L) {
        QEp.Wld <- pchisq(QE.Wld, df = QE.df, lower.tail = FALSE)
        QEp.LRT <- pchisq(QE.LRT, df = QE.df, lower.tail = FALSE)
    }
    else {
        QEp.Wld <- 1
        QEp.LRT <- 1
    }
    wi <- 1/vi
    W <- diag(wi, nrow = k.yi, ncol = k.yi)
    stXWX <- .invcalc(X = X.yi, W = W, k = k.yi)
    P <- W - W %*% X.yi %*% stXWX %*% crossprod(X.yi, W)
    vi.avg <- (k.yi - p)/.tr(P)
    I2 <- 100 * tau2/(vi.avg + tau2)
    H2 <- tau2/vi.avg + 1
    chol.h <- try(chol(vb[btt, btt]), silent = !verbose)
    if (class(chol.h) == "try-error") 
        stop("Cannot invert Hessian for QM test.")
    QM <- c(t(b)[btt] %*% chol2inv(chol.h) %*% b[btt])
    if (!int.only && int.incl && con$scale) {
        mX <- rbind(c(1, -1 * ifelse(is.d[-1], 0, meanX/sdX)), 
            cbind(0, diag(ifelse(is.d[-1], 1, 1/sdX), nrow = length(is.d) - 
                1, ncol = length(is.d) - 1)))
        b <- mX %*% b
        vb <- mX %*% vb %*% t(mX)
        X <- Xsave
    }
    rownames(vb) <- colnames(vb) <- rownames(b) <- colnames(X)
    se <- sqrt(diag(vb))
    names(se) <- NULL
    zval <- c(b/se)
    if (knha) {
        dfs <- k - p
        QM <- QM/m
        QMp <- pf(QM, df1 = m, df2 = dfs, lower.tail = FALSE)
        pval <- 2 * pt(abs(zval), df = dfs, lower.tail = FALSE)
        crit <- qt(alpha/2, df = dfs, lower.tail = FALSE)
    }
    else {
        dfs <- NA
        QMp <- pchisq(QM, df = m, lower.tail = FALSE)
        pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
        crit <- qnorm(alpha/2, lower.tail = FALSE)
    }
    ci.lb <- c(b - crit * se)
    ci.ub <- c(b + crit * se)
    if (very.verbose) 
        message("Computing fit statistics and log likelihood ...")
    ll.ML <- ifelse(method == "FE", ll.FE, ll.ML)
    ll.REML <- NA
    dev.ML <- -2 * (ll.ML - ll.QE)
    AIC.ML <- -2 * ll.ML + 2 * parms
    BIC.ML <- -2 * ll.ML + parms * log(k.eff)
    AICc.ML <- -2 * ll.ML + 2 * parms * max(k.eff, parms + 2)/(max(k.eff, 
        parms + 2) - parms - 1)
    dev.REML <- NA
    AIC.REML <- NA
    BIC.REML <- NA
    AICc.REML <- NA
    fit.stats <- matrix(c(ll.ML, dev.ML, AIC.ML, BIC.ML, AICc.ML, 
        ll.REML, dev.REML, AIC.REML, BIC.REML, AICc.REML), ncol = 2, 
        byrow = FALSE)
    dimnames(fit.stats) <- list(c("ll", "dev", "AIC", "BIC", 
        "AICc"), c("ML", "REML"))
    fit.stats <- data.frame(fit.stats)
    if (very.verbose) 
        message("Preparing output ...")
    weighted <- TRUE
    res <- list(b = b, se = se, zval = zval, pval = pval, ci.lb = ci.lb, 
        ci.ub = ci.ub, vb = vb, tau2 = tau2, se.tau2 = se.tau2, 
        sigma2 = sigma2, k = k, k.f = k.f, k.yi = k.yi, k.eff = k.eff, 
        p = p, p.eff = p.eff, parms = parms, m = m, QE.Wld = QE.Wld, 
        QEp.Wld = QEp.Wld, QE.LRT = QE.LRT, QEp.LRT = QEp.LRT, 
        QE.df = QE.df, QM = QM, QMp = QMp, I2 = I2, H2 = H2, 
        int.only = int.only, int.incl = int.incl, yi = yi, vi = vi, 
        X = X, yi.f = yi.f, vi.f = vi.f, X.f = X.f, ai = ai, 
        bi = bi, ci = ci, di = di, ai.f = ai.f, bi.f = bi.f, 
        ci.f = ci.f, di.f = di.f, x1i = x1i, x2i = x2i, t1i = t1i, 
        t2i = t2i, x1i.f = x1i.f, x2i.f = x2i.f, t1i.f = t1i.f, 
        t2i.f = t2i.f, xi = xi, mi = mi, ti = ti, xi.f = xi.f, 
        mi.f = mi.f, ti.f = ti.f, ni = ni, ni.f = ni.f, ids = ids, 
        not.na = not.na, not.na.yivi = not.na.yivi, slab = slab, 
        slab.null = slab.null, measure = measure, method = method, 
        model = model, weighted = weighted, knha = knha, dfs = dfs, 
        btt = btt, intercept = intercept, digits = digits, level = level, 
        control = control, verbose = verbose, add = add, to = to, 
        drop00 = drop00, fit.stats = fit.stats, version = packageVersion("metafor"), 
        call = mf)
    class(res) <- c("rma.glmm", "rma")
    return(res)
}
