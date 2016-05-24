rma <-
function (yi, vi, sei, weights, ai, bi, ci, di, n1i, n2i, x1i, 
    x2i, t1i, t2i, m1i, m2i, sd1i, sd2i, xi, mi, ri, ti, sdi, 
    ni, mods, measure = "GEN", intercept = TRUE, data, slab, 
    subset, add = 1/2, to = "only0", drop00 = FALSE, vtype = "LS", 
    method = "REML", weighted = TRUE, knha = FALSE, level = 95, 
    digits = 4, btt, tau2, verbose = FALSE, control) 
{
    if (!is.element(measure, c("GEN", "RR", "OR", "PETO", "RD", 
        "AS", "PHI", "YUQ", "YUY", "RTET", "PBIT", "OR2D", "OR2DN", 
        "OR2DL", "IRR", "IRD", "IRSD", "MD", "SMD", "SMDH", "ROM", 
        "RPB", "RBIS", "D2OR", "D2ORN", "D2ORL", "COR", "UCOR", 
        "ZCOR", "PR", "PLN", "PLO", "PAS", "PFT", "IR", "IRLN", 
        "IRS", "IRFT", "MN", "MC", "SMCC", "SMCR", "SMCRH", "ROMC", 
        "ARAW", "AHW", "ABT"))) 
        stop("Unknown 'measure' specified.")
    if (!is.element(method, c("FE", "HS", "HE", "DL", "GENQ", 
        "SJ", "ML", "REML", "EB", "DLIT", "SJIT", "PM"))) 
        stop("Unknown 'method' specified.")
    if (length(add) > 1) 
        add <- add[1]
    if (length(to) > 1) 
        to <- to[1]
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    if (missing(btt)) 
        btt <- NULL
    if (missing(tau2)) 
        tau2 <- NULL
    if (missing(control)) 
        control <- list()
    very.verbose <- ifelse(!is.logical(verbose) && verbose > 
        1, TRUE, FALSE)
    if (very.verbose) 
        message("Extracting/computing yi/vi values ...")
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
    mf.yi <- mf[[match("yi", names(mf))]]
    mf.weights <- mf[[match("weights", names(mf))]]
    mf.slab <- mf[[match("slab", names(mf))]]
    mf.subset <- mf[[match("subset", names(mf))]]
    mf.mods <- mf[[match("mods", names(mf))]]
    mf.scale <- mf[[match("scale", names(mf))]]
    yi <- eval(mf.yi, data, enclos = sys.frame(sys.parent()))
    weights <- eval(mf.weights, data, enclos = sys.frame(sys.parent()))
    slab <- eval(mf.slab, data, enclos = sys.frame(sys.parent()))
    subset <- eval(mf.subset, data, enclos = sys.frame(sys.parent()))
    mods <- eval(mf.mods, data, enclos = sys.frame(sys.parent()))
    scale <- eval(mf.scale, data, enclos = sys.frame(sys.parent()))
    ai <- bi <- ci <- di <- x1i <- x2i <- t1i <- t2i <- NA
    is.formula <- FALSE
    if (!is.null(yi)) {
        if (class(yi) == "formula") {
            options(na.action = "na.pass")
            mods <- model.matrix(yi, data = data)
            attr(mods, "assign") <- NULL
            yi <- model.response(model.frame(yi, data = data))
            options(na.action = na.act)
            names(yi) <- NULL
            intercept <- FALSE
            is.formula <- TRUE
        }
        if (is.matrix(yi)) 
            yi <- as.vector(yi)
        k <- length(yi)
        if (measure == "GEN" && !is.null(attr(yi, "measure"))) 
            measure <- attr(yi, "measure")
        attr(yi, "measure") <- measure
        mf.vi <- mf[[match("vi", names(mf))]]
        mf.sei <- mf[[match("sei", names(mf))]]
        mf.ni <- mf[[match("ni", names(mf))]]
        vi <- eval(mf.vi, data, enclos = sys.frame(sys.parent()))
        sei <- eval(mf.sei, data, enclos = sys.frame(sys.parent()))
        ni <- eval(mf.ni, data, enclos = sys.frame(sys.parent()))
        if (is.null(vi)) {
            if (is.null(sei)) {
                stop("Need to specify vi or sei argument.")
            }
            else {
                vi <- sei^2
            }
        }
        if (is.matrix(vi)) 
            vi <- as.vector(vi)
        if (length(vi) == 1) 
            vi <- rep(vi, k)
        if (length(vi) != k) 
            stop("Length of yi and vi (or sei) vectors are not the same.")
        if (is.null(ni) && !is.null(attr(yi, "ni"))) 
            ni <- attr(yi, "ni")
        if (!is.null(ni) && length(ni) != k) 
            ni <- NULL
        if (!is.null(ni)) 
            attr(yi, "ni") <- ni
        if (is.null(slab)) {
            if (!is.null(attr(yi, "slab"))) 
                slab <- attr(yi, "slab")
            if (is.null(slab) && length(slab) != k) 
                slab <- NULL
        }
        if (!is.null(subset)) {
            yi <- yi[subset]
            vi <- vi[subset]
            ni <- ni[subset]
            attr(yi, "measure") <- measure
            attr(yi, "ni") <- ni
        }
    }
    else {
        if (is.element(measure, c("RR", "OR", "PETO", "RD", "AS", 
            "PHI", "YUQ", "YUY", "RTET", "PBIT", "OR2D", "OR2DN", 
            "OR2DL"))) {
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
            dat <- escalc(measure = measure, ai = ai, bi = bi, 
                ci = ci, di = di, add = add, to = to, drop00 = drop00, 
                vtype = vtype)
        }
        if (is.element(measure, c("IRR", "IRD", "IRSD"))) {
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
        if (is.element(measure, c("MD", "SMD", "SMDH", "ROM", 
            "RPB", "RBIS", "D2OR", "D2ORN", "D2ORL"))) {
            mf.m1i <- mf[[match("m1i", names(mf))]]
            mf.m2i <- mf[[match("m2i", names(mf))]]
            mf.sd1i <- mf[[match("sd1i", names(mf))]]
            mf.sd2i <- mf[[match("sd2i", names(mf))]]
            mf.n1i <- mf[[match("n1i", names(mf))]]
            mf.n2i <- mf[[match("n2i", names(mf))]]
            m1i <- eval(mf.m1i, data, enclos = sys.frame(sys.parent()))
            m2i <- eval(mf.m2i, data, enclos = sys.frame(sys.parent()))
            sd1i <- eval(mf.sd1i, data, enclos = sys.frame(sys.parent()))
            sd2i <- eval(mf.sd2i, data, enclos = sys.frame(sys.parent()))
            n1i <- eval(mf.n1i, data, enclos = sys.frame(sys.parent()))
            n2i <- eval(mf.n2i, data, enclos = sys.frame(sys.parent()))
            k <- length(m1i)
            if (!is.null(subset)) {
                m1i <- m1i[subset]
                m2i <- m2i[subset]
                sd1i <- sd1i[subset]
                sd2i <- sd2i[subset]
                n1i <- n1i[subset]
                n2i <- n2i[subset]
            }
            dat <- escalc(measure = measure, m1i = m1i, m2i = m2i, 
                sd1i = sd1i, sd2i = sd2i, n1i = n1i, n2i = n2i, 
                vtype = vtype)
        }
        if (is.element(measure, c("COR", "UCOR", "ZCOR"))) {
            mf.ri <- mf[[match("ri", names(mf))]]
            mf.ni <- mf[[match("ni", names(mf))]]
            ri <- eval(mf.ri, data, enclos = sys.frame(sys.parent()))
            ni <- eval(mf.ni, data, enclos = sys.frame(sys.parent()))
            k <- length(ri)
            if (!is.null(subset)) {
                ri <- ri[subset]
                ni <- ni[subset]
            }
            dat <- escalc(measure = measure, ri = ri, ni = ni, 
                vtype = vtype)
        }
        if (is.element(measure, c("PR", "PLN", "PLO", "PAS", 
            "PFT"))) {
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
            dat <- escalc(measure = measure, xi = xi, mi = mi, 
                add = add, to = to, vtype = vtype)
        }
        if (is.element(measure, c("IR", "IRLN", "IRS", "IRFT"))) {
            mf.xi <- mf[[match("xi", names(mf))]]
            mf.ti <- mf[[match("ti", names(mf))]]
            xi <- eval(mf.xi, data, enclos = sys.frame(sys.parent()))
            ti <- eval(mf.ti, data, enclos = sys.frame(sys.parent()))
            k <- length(xi)
            if (!is.null(subset)) {
                xi <- xi[subset]
                ti <- ti[subset]
            }
            dat <- escalc(measure = measure, xi = xi, ti = ti, 
                add = add, to = to, vtype = vtype)
        }
        if (is.element(measure, c("MN"))) {
            mf.mi <- mf[[match("mi", names(mf))]]
            mf.sdi <- mf[[match("sdi", names(mf))]]
            mf.ni <- mf[[match("ni", names(mf))]]
            mi <- eval(mf.mi, data, enclos = sys.frame(sys.parent()))
            sdi <- eval(mf.sdi, data, enclos = sys.frame(sys.parent()))
            ni <- eval(mf.ni, data, enclos = sys.frame(sys.parent()))
            k <- length(mi)
            if (!is.null(subset)) {
                mi <- mi[subset]
                sdi <- sdi[subset]
                ni <- ni[subset]
            }
            dat <- escalc(measure = measure, mi = mi, sdi = sdi, 
                ni = ni, vtype = vtype)
        }
        if (is.element(measure, c("MC", "SMCC", "SMCR", "SMCRH", 
            "ROMC"))) {
            mf.m1i <- mf[[match("m1i", names(mf))]]
            mf.m2i <- mf[[match("m2i", names(mf))]]
            mf.sd1i <- mf[[match("sd1i", names(mf))]]
            mf.sd2i <- mf[[match("sd2i", names(mf))]]
            mf.ri <- mf[[match("ri", names(mf))]]
            mf.ni <- mf[[match("ni", names(mf))]]
            m1i <- eval(mf.m1i, data, enclos = sys.frame(sys.parent()))
            m2i <- eval(mf.m2i, data, enclos = sys.frame(sys.parent()))
            sd1i <- eval(mf.sd1i, data, enclos = sys.frame(sys.parent()))
            sd2i <- eval(mf.sd2i, data, enclos = sys.frame(sys.parent()))
            ri <- eval(mf.ri, data, enclos = sys.frame(sys.parent()))
            ni <- eval(mf.ni, data, enclos = sys.frame(sys.parent()))
            k <- length(m1i)
            if (!is.null(subset)) {
                m1i <- m1i[subset]
                m2i <- m2i[subset]
                sd1i <- sd1i[subset]
                sd2i <- sd2i[subset]
                ni <- ni[subset]
                ri <- ri[subset]
            }
            dat <- escalc(measure = measure, m1i = m1i, m2i = m2i, 
                sd1i = sd1i, sd2i = sd2i, ri = ri, ni = ni, vtype = vtype)
        }
        if (is.element(measure, c("ARAW", "AHW", "ABT"))) {
            mf.ai <- mf[[match("ai", names(mf))]]
            mf.mi <- mf[[match("mi", names(mf))]]
            mf.ni <- mf[[match("ni", names(mf))]]
            ai <- eval(mf.ai, data, enclos = sys.frame(sys.parent()))
            mi <- eval(mf.mi, data, enclos = sys.frame(sys.parent()))
            ni <- eval(mf.ni, data, enclos = sys.frame(sys.parent()))
            k <- length(ai)
            if (!is.null(subset)) {
                ai <- ai[subset]
                mi <- mi[subset]
                ni <- ni[subset]
            }
            dat <- escalc(measure = measure, ai = ai, mi = mi, 
                ni = ni, vtype = vtype)
        }
        if (is.element(measure, "GEN")) 
            stop("Specify the desired outcome measure via the 'measure' argument.")
        yi <- dat$yi
        vi <- dat$vi
        ni <- attr(yi, "ni")
    }
    if (length(weights) == 1) 
        weights <- rep(weights, k)
    if (!is.null(weights) && (length(weights) != k)) 
        stop("Length of yi and weights vectors are not the same.")
    if (!is.null(subset)) 
        weights <- weights[subset]
    if (very.verbose) 
        message("Creating model matrix ...")
    if (class(mods) == "formula") {
        options(na.action = "na.pass")
        mods <- model.matrix(mods, data = data)
        attr(mods, "assign") <- NULL
        options(na.action = na.act)
        intercept <- FALSE
        is.formula <- TRUE
    }
    if (is.vector(mods)) 
        mods <- cbind(mods)
    if (is.data.frame(mods)) 
        mods <- as.matrix(mods)
    if (is.character(mods)) 
        stop("Model matrix contains character variables.")
    if (!is.null(mods) && (nrow(mods) != k)) 
        stop("Number of rows of the model matrix does not match length of yi argument.")
    if (class(scale) == "formula") {
        options(na.action = "na.pass")
        Z <- model.matrix(scale, data = data)
        attr(Z, "assign") <- NULL
        options(na.action = na.act)
        model <- "rma.tau2"
        if (nrow(Z) != k) 
            stop("Number of rows of the model matrix for tau2 does not match length of yi argument.")
    }
    else {
        Z <- NULL
        model <- "rma.uni"
    }
    ids <- seq_len(k)
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
        Z <- Z[subset, , drop = FALSE]
    }
    if (anyDuplicated(slab)) 
        slab <- make.unique(as.character(slab))
    attr(yi, "slab") <- slab
    k <- length(yi)
    if (any(vi <= 0, na.rm = TRUE)) {
        allvipos <- FALSE
        warning("There are outcomes with non-positive sampling variances.")
        vi.neg <- vi < 0
        if (any(vi.neg, na.rm = TRUE)) {
            vi[vi.neg] <- 0
            warning("Negative sampling variances constrained to zero.")
        }
    }
    else {
        allvipos <- TRUE
    }
    if (any(weights < 0, na.rm = TRUE)) 
        stop("Negative weights not allowed.")
    if (any(is.infinite(weights))) 
        stop("Infinite weights not allowed.")
    ai.f <- ai
    bi.f <- bi
    ci.f <- ci
    di.f <- di
    x1i.f <- x1i
    x2i.f <- x2i
    t1i.f <- t1i
    t2i.f <- t2i
    yi.f <- yi
    vi.f <- vi
    weights.f <- weights
    ni.f <- ni
    mods.f <- mods
    k.f <- k
    YVXZW.na <- is.na(yi) | is.na(vi) | if (is.null(mods)) 
        FALSE
    else apply(is.na(mods), 1, any) | if (is.null(Z)) 
        FALSE
    else apply(is.na(Z), 1, any) | if (is.null(weights)) 
        FALSE
    else is.na(weights)
    if (any(YVXZW.na)) {
        if (very.verbose) 
            message("Handling NAs ...")
        not.na <- !YVXZW.na
        if (na.act == "na.omit" || na.act == "na.exclude" || 
            na.act == "na.pass") {
            yi <- yi[not.na]
            vi <- vi[not.na]
            weights <- weights[not.na]
            ni <- ni[not.na]
            mods <- mods[not.na, , drop = FALSE]
            Z <- Z[not.na, , drop = FALSE]
            k <- length(yi)
            warning("Studies with NAs omitted from model fitting.")
            attr(yi, "measure") <- measure
            attr(yi, "ni") <- ni
        }
        if (na.act == "na.fail") 
            stop("Missing values in data.")
    }
    else {
        not.na <- rep(TRUE, k)
    }
    if (k < 1) 
        stop("Processing terminated since k = 0.")
    if (k == 1) {
        method <- "FE"
        knha <- FALSE
    }
    if (is.null(mods) && !intercept) {
        warning("Must either include an intercept and/or moderators in model.\n  Coerced intercept into the model.")
        intercept <- TRUE
    }
    if (intercept) {
        X <- cbind(intrcpt = rep(1, k), mods)
        X.f <- cbind(intrcpt = rep(1, k.f), mods.f)
    }
    else {
        X <- mods
        X.f <- mods.f
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
    tmp <- lm(yi ~ X - 1)
    coef.na <- is.na(coef(tmp))
    if (any(coef.na)) {
        warning("Redundant predictors dropped from the model.")
        X <- X[, !coef.na, drop = FALSE]
        X.f <- X.f[, !coef.na, drop = FALSE]
    }
    p <- NCOL(X)
    if ((p == 1L) && (all(sapply(X, identical, 1)))) {
        int.only <- TRUE
    }
    else {
        int.only <- FALSE
    }
    if (method == "FE") {
        if (p > k) 
            stop("Number of parameters to be estimated is larger than the number of observations.")
    }
    else {
        if (is.numeric(tau2)) {
            if (p > k) 
                stop("Number of parameters to be estimated is larger than the number of observations.")
        }
        else {
            if ((p + 1) > k) 
                stop("Number of parameters to be estimated is larger than the number of observations.")
        }
    }
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
    con <- list(verbose = FALSE, tau2.init = NULL, tau2.min = 0, 
        tau2.max = 100, threshold = 10^-5, maxiter = 100, stepadj = 1, 
        REMLf = TRUE, tol = 1e-07, optimizer = "nlminb", optmethod = "BFGS")
    con.pos <- pmatch(names(control), names(con))
    con[c(na.omit(con.pos))] <- control[!is.na(con.pos)]
    if (verbose) 
        con$verbose <- verbose
    verbose <- con$verbose
    if (any(eigen(crossprod(X), symmetric = TRUE, only.values = TRUE)$values <= 
        con$tol)) 
        stop("Model matrix not of full rank. Cannot fit model.")
    iter <- 0
    se.tau2 <- I2 <- H2 <- QE <- QEp <- NA
    s2w <- 1
    alpha <- ifelse(level > 1, (100 - level)/100, 1 - level)
    Y <- as.matrix(yi)
    if (model == "rma.tau2") {
        if (method != "ML" && method != "REML") 
            stop("Location-scale model can only be fitted with ML or REML estimation.")
        if (!weighted) 
            stop("Cannot use weighted estimation to fit location-scale model.")
        if (!is.null(weights)) 
            stop("Cannot use user-defined weights for location-scale model.")
        if (any(eigen(crossprod(Z), symmetric = TRUE, only.values = TRUE)$values <= 
            con$tol)) 
            stop("Model matrix for scale part of the model not of full rank. Cannot fit model.")
        optimizer <- match.arg(con$optimizer, c("optim", "nlminb", 
            "uobyqa", "newuoa", "bobyqa", "nloptr"))
        optmethod <- con$optmethod
        optcontrol <- control[is.na(con.pos)]
        if (length(optcontrol) == 0) 
            optcontrol <- list()
        if (optimizer == "nloptr" && !is.element("algorithm", 
            names(optcontrol))) 
            optcontrol$algorithm <- "NLOPT_LN_BOBYQA"
        if (optimizer == "nloptr" && !is.element("ftol_rel", 
            names(optcontrol))) 
            optcontrol$ftol_rel <- 1e-08
        reml <- ifelse(method == "REML", TRUE, FALSE)
        if (is.element(optimizer, c("uobyqa", "newuoa", "bobyqa"))) {
            if (!requireNamespace("minqa", quietly = TRUE)) 
                stop("Please install the 'minqa' package to use this optimizer.")
        }
        if (optimizer == "nloptr") {
            if (!requireNamespace("nloptr", quietly = TRUE)) 
                stop("Please install the 'nloptr' package to use this optimizer.")
        }
        if (!requireNamespace("numDeriv", quietly = TRUE)) 
            stop("Please install the 'numDeriv' package to fit a location-scale model.")
        p.tau2 <- NCOL(Z)
        if (!is.null(con$tau2.init)) {
            if (length(con$tau2.init) != p.tau2) 
                stop(paste("Length of 'tau2.init' argument (", 
                  length(con$tau2.init), ") does not match actual number of parameters (", 
                  p.tau2, ").", sep = ""))
        }
        else {
            con$tau2.init <- rep(1e-04, p.tau2)
        }
        if (optimizer == "optim") 
            par.arg <- "par"
        if (optimizer == "nlminb") 
            par.arg <- "start"
        if (is.element(optimizer, c("uobyqa", "newuoa", "bobyqa"))) {
            par.arg <- "par"
            optimizer <- paste0("minqa::", optimizer)
        }
        if (optimizer == "nloptr") {
            par.arg <- "x0"
            optimizer <- paste0("nloptr::nloptr")
        }
        optcall <- paste(optimizer, "(", par.arg, "=con$tau2.init, .ll.rma.tau2, ", 
            ifelse(optimizer == "optim", "method=optmethod, ", 
                ""), "yi=yi, vi=vi, X=X, Z=Z, reml=reml,\n                                                    k=k, p=p,\n                                                    verbose=verbose, digits=digits, REMLf=con$REMLf, ", 
            ifelse(optimizer == "nloptr::nloptr", "opts=optcontrol)", 
                "control=optcontrol)"), sep = "")
        opt.res <- try(eval(parse(text = optcall)), silent = !verbose)
        if (inherits(opt.res, "try-error")) 
            stop("Error during optimization for location-scale model.")
        if (is.element(optimizer, c("optim", "nlminb")) && opt.res$convergence != 
            0) 
            stop(paste0("Optimizer (", optimizer, ") did not achieve convergence (convergence = ", 
                opt.res$convergence, ")."))
        if (is.element(optimizer, c("minqa::uobyqa", "minqa::newuoa", 
            "minqa::bobyqa")) && opt.res$ierr != 0) 
            stop(paste0("Optimizer (", optimizer, ") did not achieve convergence (ierr = ", 
                opt.res$ierr, ")."))
        if (optimizer == "nloptr::nloptr" && !(opt.res$status >= 
            1 && opt.res$status <= 4)) 
            stop(paste0("Optimizer (", optimizer, ") did not achieve convergence (status = ", 
                opt.res$status, ")."))
        if (verbose) {
            cat("\n")
            print(opt.res)
        }
        if (optimizer == "nloptr::nloptr") 
            opt.res$par <- opt.res$solution
        vb.tau2 <- matrix(NA, nrow = p.tau2, ncol = p.tau2)
        se.tau2 <- rep(NA, p.tau2)
        h <- try(numDeriv::hessian(.ll.rma.tau2, x = opt.res$par, 
            yi = yi, vi = vi, X = X, Z = Z, reml = reml, k = k, 
            p = p, verbose = FALSE, digits = digits, REMLf = con$REMLf), 
            silent = TRUE)
        if (!inherits(h, "try-error")) {
            chol.h <- try(chol(h), silent = TRUE)
            if (!inherits(chol.h, "try-error")) {
                vb.tau2 <- chol2inv(chol.h)
                se.tau2 <- sqrt(diag(vb.tau2))
            }
        }
        b.tau2 <- cbind(opt.res$par)
        colnames(Z)[grep("(Intercept)", colnames(Z))] <- "intrcpt"
        rownames(b.tau2) <- rownames(vb.tau2) <- colnames(vb.tau2) <- colnames(Z)
        crit <- qnorm(alpha/2, lower.tail = FALSE)
        names(se.tau2) <- NULL
        zval.tau2 <- c(b.tau2/se.tau2)
        pval.tau2 <- 2 * pnorm(abs(zval.tau2), lower.tail = FALSE)
        ci.lb.tau2 <- c(b.tau2 - crit * se.tau2)
        ci.ub.tau2 <- c(b.tau2 + crit * se.tau2)
        tau2 <- exp(as.vector(Z %*% b.tau2))
        tau2.fix <- FALSE
    }
    else {
        if (is.numeric(tau2)) {
            tau2.fix <- TRUE
            tau2.val <- tau2
        }
        else {
            tau2.fix <- FALSE
            tau2.val <- NA
        }
        if (very.verbose && !tau2.fix) 
            message("Estimating tau^2 value ...")
        if (method == "HS") {
            if (!allvipos) 
                stop("HS estimator cannot be used with non-positive sampling variances.")
            wi <- 1/vi
            W <- diag(wi, nrow = k, ncol = k)
            stXWX <- .invcalc(X = X, W = W, k = k)
            P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
            RSS <- crossprod(Y, P) %*% Y
            tau2 <- ifelse(tau2.fix, tau2.val, (RSS - k)/sum(wi))
        }
        if (is.element(method, c("HE", "ML", "REML", "EB"))) {
            stXX <- .invcalc(X = X, W = diag(k), k = k)
            P <- diag(k) - X %*% tcrossprod(stXX, X)
            RSS <- crossprod(Y, P) %*% Y
            V <- diag(vi, nrow = k, ncol = k)
            PV <- P %*% V
            trPV <- .tr(PV)
            tau2 <- ifelse(tau2.fix, tau2.val, (RSS - trPV)/(k - 
                p))
        }
        if (method == "DL") {
            if (!allvipos) 
                stop("DL estimator cannot be used with non-positive sampling variances.")
            wi <- 1/vi
            W <- diag(wi, nrow = k, ncol = k)
            stXWX <- .invcalc(X = X, W = W, k = k)
            P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
            RSS <- crossprod(Y, P) %*% Y
            trP <- .tr(P)
            tau2 <- ifelse(tau2.fix, tau2.val, (RSS - (k - p))/trP)
        }
        if (method == "GENQ") {
            if (is.null(weights)) 
                stop("Must specify 'weights' when method='GENQ'.")
            A <- diag(weights, nrow = k, ncol = k)
            stXAX <- .invcalc(X = X, W = A, k = k)
            P <- A - A %*% X %*% stXAX %*% t(X) %*% A
            RSS <- crossprod(Y, P) %*% Y
            V <- diag(vi, nrow = k, ncol = k)
            PV <- P %*% V
            trP <- .tr(P)
            trPV <- .tr(PV)
            tau2 <- ifelse(tau2.fix, tau2.val, (RSS - trPV)/trP)
        }
        if (method == "SJ") {
            if (is.null(con$tau2.init)) {
                tau2.0 <- c(var(yi) * (k - 1)/k)
            }
            else {
                tau2.0 <- con$tau2.init
            }
            wi <- 1/(vi + tau2.0)
            W <- diag(wi, nrow = k, ncol = k)
            stXWX <- .invcalc(X = X, W = W, k = k)
            P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
            RSS <- crossprod(Y, P) %*% Y
            V <- diag(vi, nrow = k, ncol = k)
            PV <- P %*% V
            tau2 <- ifelse(tau2.fix, tau2.val, tau2.0 * RSS/(k - 
                p))
        }
        if (method == "DLIT") {
            conv <- 1
            change <- con$threshold + 1
            if (is.null(con$tau2.init)) {
                tau2 <- 0
            }
            else {
                tau2 <- con$tau2.init
            }
            while (change > con$threshold) {
                if (verbose) 
                  cat("Iteration", iter, "\ttau^2 =", formatC(tau2, 
                    format = "f", digits = digits), "\n")
                iter <- iter + 1
                tau2.old <- tau2
                wi <- 1/(vi + tau2)
                if (any(is.infinite(wi))) 
                  stop("Division by zero when computing the inverse variance weights.")
                W <- diag(wi, nrow = k, ncol = k)
                stXWX <- .invcalc(X = X, W = W, k = k)
                P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
                RSS <- crossprod(Y, P) %*% Y
                trP <- .tr(P)
                tau2 <- ifelse(tau2.fix, tau2.val, (RSS - (k - 
                  p))/trP)
                tau2[tau2 < con$tau2.min] <- con$tau2.min
                change <- abs(tau2.old - tau2)
                if (iter > con$maxiter) {
                  conv <- 0
                  break
                }
            }
            if (conv == 0L) 
                stop("Algorithm did not converge.")
        }
        if (method == "SJIT") {
            conv <- 1
            change <- con$threshold + 1
            if (is.null(con$tau2.init)) {
                tau2 <- var(yi) * (k - 1)/k
                tau2.0 <- tau2
            }
            else {
                tau2 <- con$tau2.init
                tau2.0 <- tau2
            }
            while (change > con$threshold) {
                if (verbose) 
                  cat("Iteration", iter, "\ttau^2 =", formatC(tau2, 
                    format = "f", digits = digits), "\n")
                iter <- iter + 1
                tau2.old <- tau2
                wi <- 1/(vi + tau2)
                W <- diag(wi, nrow = k, ncol = k)
                stXWX <- .invcalc(X = X, W = W, k = k)
                P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
                RSS <- crossprod(Y, P) %*% Y
                V <- diag(vi, nrow = k, ncol = k)
                PV <- P %*% V
                tau2 <- ifelse(tau2.fix, tau2.val, tau2 * RSS/(k - 
                  p))
                change <- abs(tau2.old - tau2)
                if (iter > con$maxiter) {
                  conv <- 0
                  break
                }
            }
            if (conv == 0L) 
                stop("Algorithm did not converge.")
        }
        if (method == "PM") {
            if (!allvipos) 
                stop("PM estimator cannot be used with non-positive sampling variances.")
            if (.QE.func(con$tau2.min, Y = Y, vi = vi, X = X, 
                k = k, objective = k - p) < con$tau2.min) {
                tau2 <- con$tau2.min
            }
            else {
                tau2 <- ifelse(tau2.fix, tau2.val, try(uniroot(.QE.func, 
                  interval = c(con$tau2.min, con$tau2.max), tol = con$threshold, 
                  maxiter = con$maxiter, Y = Y, vi = vi, X = X, 
                  k = k, objective = k - p, verbose = verbose, 
                  digits = digits, extendInt = "upX")$root, silent = TRUE))
                if (!is.numeric(tau2)) 
                  stop("Error in iterative search for tau2. Try increasing tau2.max or switch to another 'method'.")
            }
            wi <- 1/(vi + tau2)
            W <- diag(wi, nrow = k, ncol = k)
            stXWX <- .invcalc(X = X, W = W, k = k)
            P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
        }
        if (is.element(method, c("ML", "REML", "EB"))) {
            conv <- 1
            change <- con$threshold + 1
            if (is.null(con$tau2.init)) {
                tau2 <- max(0, tau2, con$tau2.min)
            }
            else {
                tau2 <- con$tau2.init
            }
            while (change > con$threshold) {
                if (verbose) 
                  cat("Iteration", iter, "\ttau^2 =", formatC(tau2, 
                    format = "f", digits = digits), "\n")
                iter <- iter + 1
                tau2.old <- tau2
                wi <- 1/(vi + tau2)
                if (any(is.infinite(wi))) 
                  stop("Division by zero when computing the inverse variance weights.")
                W <- diag(wi, nrow = k, ncol = k)
                stXWX <- .invcalc(X = X, W = W, k = k)
                P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
                if (method == "ML") {
                  PP <- P %*% P
                  adj <- (crossprod(Y, PP) %*% Y - sum(wi))/sum(wi^2)
                }
                if (method == "REML") {
                  PP <- P %*% P
                  adj <- (crossprod(Y, PP) %*% Y - .tr(P))/.tr(PP)
                }
                if (method == "EB") {
                  adj <- (crossprod(Y, P) %*% Y * k/(k - p) - 
                    k)/sum(wi)
                }
                adj <- adj * con$stepadj
                while (tau2 + adj < con$tau2.min) {
                  adj <- adj/2
                }
                tau2 <- ifelse(tau2.fix, tau2.val, tau2 + adj)
                change <- abs(tau2.old - tau2)
                if (iter > con$maxiter) {
                  conv <- 0
                  break
                }
            }
            if (conv == 0L) 
                stop("Fisher scoring algorithm did not converge. See 'help(rma)' for possible remedies.")
        }
        tau2 <- max(con$tau2.min, c(tau2))
        if (verbose && is.element(method, c("ML", "REML", "EB"))) {
            cat("Iteration", iter, "\ttau^2 =", formatC(tau2, 
                format = "f", digits = digits), "\n")
            cat("Fisher scoring algorithm converged after", iter, 
                "iterations.\n")
        }
        if (method == "HS") {
            se.tau2 <- sqrt(1/sum(wi)^2 * (2 * (k - p) + 4 * 
                max(tau2, 0) * .tr(P) + 2 * max(tau2, 0)^2 * 
                sum(P * P)))
        }
        if (method == "HE") {
            se.tau2 <- sqrt(1/(k - p)^2 * (2 * sum(PV * t(PV)) + 
                4 * max(tau2, 0) * trPV + 2 * max(tau2, 0)^2 * 
                (k - p)))
        }
        if (method == "DL" || method == "DLIT") {
            se.tau2 <- sqrt(1/trP^2 * (2 * (k - p) + 4 * max(tau2, 
                0) * trP + 2 * max(tau2, 0)^2 * sum(P * P)))
        }
        if (method == "GENQ") {
            se.tau2 <- sqrt(1/trP^2 * (2 * sum(PV * t(PV)) + 
                4 * max(tau2, 0) * sum(PV * P) + 2 * max(tau2, 
                0)^2 * sum(P * P)))
        }
        if (method == "SJ") {
            se.tau2 <- sqrt(tau2.0^2/(k - p)^2 * (2 * sum(PV * 
                t(PV)) + 4 * max(tau2, 0) * sum(PV * P) + 2 * 
                max(tau2, 0)^2 * sum(P * P)))
        }
        if (method == "ML") {
            se.tau2 <- sqrt(2/sum(wi^2))
        }
        if (method == "REML") {
            se.tau2 <- sqrt(2/sum(P * P))
        }
        if (method == "EB" || method == "PM" || method == "SJIT") {
            V <- diag(vi, nrow = k, ncol = k)
            PV <- P %*% V
            se.tau2 <- sqrt(2 * k^2/(k - p)/sum(wi)^2)
        }
    }
    if (method == "FE") 
        tau2 <- 0
    if (very.verbose) 
        message("Model fitting ...")
    wi <- 1/(vi + tau2)
    W <- diag(wi, nrow = k, ncol = k)
    M <- diag(vi + tau2, nrow = k, ncol = k)
    if (weighted) {
        if (is.null(weights)) {
            if (any(is.infinite(wi))) 
                stop("Division by zero when computing the inverse variance weights.")
            stXWX <- .invcalc(X = X, W = W, k = k)
            b <- stXWX %*% crossprod(X, W) %*% Y
            vb <- stXWX
            RSS.f <- sum(wi * (yi - X %*% b)^2)
        }
        else {
            A <- diag(weights, nrow = k, ncol = k)
            stXAX <- .invcalc(X = X, W = A, k = k)
            b <- stXAX %*% crossprod(X, A) %*% Y
            vb <- stXAX %*% t(X) %*% A %*% M %*% A %*% X %*% 
                stXAX
            RSS.f <- sum(wi * (yi - X %*% b)^2)
        }
        if (knha) {
            if (RSS.f <= .Machine$double.eps) {
                s2w <- 1
            }
            else {
                s2w <- RSS.f/(k - p)
            }
            vb <- s2w * vb
            if (method == "FE") 
                warning("Knapp & Hartung (2003) method is not meant to be used in the context of fixed-effects models.")
        }
        QM <- c(t(b)[btt] %*% chol2inv(chol(vb[btt, btt])) %*% 
            b[btt])
    }
    else {
        stXX <- .invcalc(X = X, W = diag(k), k = k)
        b <- stXX %*% crossprod(X, Y)
        vb <- tcrossprod(stXX, X) %*% M %*% X %*% stXX
        RSS.f <- sum(wi * (yi - X %*% b)^2)
        if (knha) {
            stXWX <- .invcalc(X = X, W = W, k = k)
            b.knha <- stXWX %*% crossprod(X, W) %*% Y
            RSS.knha <- sum(wi * (yi - X %*% b.knha)^2)
            if (RSS.knha <= .Machine$double.eps) {
                s2w <- 1
            }
            else {
                s2w <- RSS.knha/(k - p)
            }
            vb <- s2w * vb
            if (method == "FE") 
                warning("Knapp & Hartung (2003) method is not meant to be used in the context of fixed-effects models.")
        }
        QM <- c(t(b)[btt] %*% chol2inv(chol(vb[btt, btt])) %*% 
            b[btt])
    }
    rownames(b) <- rownames(vb) <- colnames(vb) <- colnames(X)
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
        message("Heterogeneity testing ...")
    if (allvipos) {
        wi <- 1/vi
        W.FE <- diag(wi, nrow = k, ncol = k)
        stXWX <- .invcalc(X = X, W = W.FE, k = k)
        P <- W.FE - W.FE %*% X %*% stXWX %*% crossprod(X, W.FE)
        QE <- max(0, c(crossprod(Y, P) %*% Y))
        QEp <- ifelse(k - p >= 1, pchisq(QE, df = k - p, lower.tail = FALSE), 
            1)
        vi.avg <- (k - p)/.tr(P)
        I2 <- 100 * mean(tau2)/(vi.avg + mean(tau2))
        H2 <- mean(tau2)/vi.avg + 1
    }
    else {
        warning(paste0("Cannot compute ", ifelse(int.only, "Q", 
            "QE"), "-test, I^2, or H^2 with non-positive sampling variances."))
    }
    if (!int.only && int.incl && method != "FE" && model == "rma.uni") {
        if (very.verbose) {
            message("Fitting RE model for R^2 computation ...")
            res.RE <- try(rma(yi, vi, weights = weights, method = method, 
                weighted = weighted, knha = knha, verbose = ifelse(verbose, 
                  TRUE, FALSE), control = con, digits = digits), 
                silent = FALSE)
        }
        else {
            res.RE <- suppressWarnings(try(rma(yi, vi, weights = weights, 
                method = method, weighted = weighted, knha = knha, 
                verbose = ifelse(verbose, TRUE, FALSE), control = con, 
                digits = digits), silent = FALSE))
        }
        if (!inherits(res.RE, "try-error")) {
            tau2.RE <- res.RE$tau2
            if (identical(tau2.RE, 0)) {
                R2 <- NA
            }
            else {
                R2 <- round(max(0, 100 * (tau2.RE - tau2)/tau2.RE), 
                  2)
            }
        }
        else {
            R2 <- NA
        }
    }
    else {
        R2 <- NULL
    }
    if (very.verbose) 
        message("Computing fit statistics and log likelihood ...")
    parms <- p + ifelse(model == "rma.uni", ifelse(method == 
        "FE" || tau2.fix, 0, 1), p.tau2)
    ll.ML <- -1/2 * (k) * log(2 * base::pi) - 1/2 * sum(log(vi + 
        tau2)) - 1/2 * RSS.f
    ll.REML <- -1/2 * (k - p) * log(2 * base::pi) + ifelse(con$REMLf, 
        1/2 * determinant(crossprod(X), logarithm = TRUE)$modulus, 
        0) - 1/2 * sum(log(vi + tau2)) - 1/2 * determinant(crossprod(X, 
        W) %*% X, logarithm = TRUE)$modulus - 1/2 * RSS.f
    if (k - p > 0L) {
        dev.ML <- -2 * (ll.ML - sum(dnorm(yi, mean = yi, sd = sqrt(vi), 
            log = TRUE)))
    }
    else {
        dev.ML <- 0
    }
    AIC.ML <- -2 * ll.ML + 2 * parms
    BIC.ML <- -2 * ll.ML + parms * log(k)
    AICc.ML <- -2 * ll.ML + 2 * parms * max(k, parms + 2)/(max(k, 
        parms + 2) - parms - 1)
    dev.REML <- -2 * (ll.REML - 0)
    AIC.REML <- -2 * ll.REML + 2 * parms
    BIC.REML <- -2 * ll.REML + parms * log(k - p)
    AICc.REML <- -2 * ll.REML + 2 * parms * max(k - p, parms + 
        2)/(max(k - p, parms + 2) - parms - 1)
    fit.stats <- matrix(c(ll.ML, dev.ML, AIC.ML, BIC.ML, AICc.ML, 
        ll.REML, dev.REML, AIC.REML, BIC.REML, AICc.REML), ncol = 2, 
        byrow = FALSE)
    dimnames(fit.stats) <- list(c("ll", "dev", "AIC", "BIC", 
        "AICc"), c("ML", "REML"))
    fit.stats <- data.frame(fit.stats)
    if (very.verbose) 
        message("Preparing output ...")
    p.eff <- p
    k.eff <- k
    res <- list(b = b, se = se, zval = zval, pval = pval, ci.lb = ci.lb, 
        ci.ub = ci.ub, vb = vb, tau2 = tau2, se.tau2 = se.tau2, 
        tau2.fix = tau2.fix, k = k, k.f = k.f, k.eff = k.eff, 
        p = p, p.eff = p.eff, parms = parms, m = m, QE = QE, 
        QEp = QEp, QM = QM, QMp = QMp, I2 = I2, H2 = H2, R2 = R2, 
        int.only = int.only, int.incl = int.incl, allvipos = allvipos, 
        coef.na = coef.na, yi = yi, vi = vi, X = X, weights = weights, 
        yi.f = yi.f, vi.f = vi.f, X.f = X.f, weights.f = weights.f, 
        M = M, ai.f = ai.f, bi.f = bi.f, ci.f = ci.f, di.f = di.f, 
        x1i.f = x1i.f, x2i.f = x2i.f, t1i.f = t1i.f, t2i.f = t2i.f, 
        ni = ni, ni.f = ni.f, ids = ids, not.na = not.na, subset = subset, 
        slab = slab, slab.null = slab.null, measure = measure, 
        method = method, weighted = weighted, knha = knha, dfs = dfs, 
        s2w = s2w, btt = btt, intercept = intercept, digits = digits, 
        level = level, sparse = FALSE, control = control, verbose = verbose, 
        add = add, to = to, drop00 = drop00, fit.stats = fit.stats, 
        version = packageVersion("metafor"), model = model, call = mf)
    if (model == "rma.tau2") {
        res$b.tau2 <- b.tau2
        res$vb.tau2 <- vb.tau2
        res$se.tau2 <- se.tau2
        res$zval.tau2 <- zval.tau2
        res$pval.tau2 <- pval.tau2
        res$ci.lb.tau2 <- ci.lb.tau2
        res$ci.ub.tau2 <- ci.ub.tau2
        res$Z <- Z
    }
    class(res) <- c("rma.uni", "rma")
    return(res)
}
rma.uni <-
function (yi, vi, sei, weights, ai, bi, ci, di, n1i, n2i, x1i, 
    x2i, t1i, t2i, m1i, m2i, sd1i, sd2i, xi, mi, ri, ti, sdi, 
    ni, mods, measure = "GEN", intercept = TRUE, data, slab, 
    subset, add = 1/2, to = "only0", drop00 = FALSE, vtype = "LS", 
    method = "REML", weighted = TRUE, knha = FALSE, level = 95, 
    digits = 4, btt, tau2, verbose = FALSE, control) 
{
    if (!is.element(measure, c("GEN", "RR", "OR", "PETO", "RD", 
        "AS", "PHI", "YUQ", "YUY", "RTET", "PBIT", "OR2D", "OR2DN", 
        "OR2DL", "IRR", "IRD", "IRSD", "MD", "SMD", "SMDH", "ROM", 
        "RPB", "RBIS", "D2OR", "D2ORN", "D2ORL", "COR", "UCOR", 
        "ZCOR", "PR", "PLN", "PLO", "PAS", "PFT", "IR", "IRLN", 
        "IRS", "IRFT", "MN", "MC", "SMCC", "SMCR", "SMCRH", "ROMC", 
        "ARAW", "AHW", "ABT"))) 
        stop("Unknown 'measure' specified.")
    if (!is.element(method, c("FE", "HS", "HE", "DL", "GENQ", 
        "SJ", "ML", "REML", "EB", "DLIT", "SJIT", "PM"))) 
        stop("Unknown 'method' specified.")
    if (length(add) > 1) 
        add <- add[1]
    if (length(to) > 1) 
        to <- to[1]
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    if (missing(btt)) 
        btt <- NULL
    if (missing(tau2)) 
        tau2 <- NULL
    if (missing(control)) 
        control <- list()
    very.verbose <- ifelse(!is.logical(verbose) && verbose > 
        1, TRUE, FALSE)
    if (very.verbose) 
        message("Extracting/computing yi/vi values ...")
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
    mf.yi <- mf[[match("yi", names(mf))]]
    mf.weights <- mf[[match("weights", names(mf))]]
    mf.slab <- mf[[match("slab", names(mf))]]
    mf.subset <- mf[[match("subset", names(mf))]]
    mf.mods <- mf[[match("mods", names(mf))]]
    mf.scale <- mf[[match("scale", names(mf))]]
    yi <- eval(mf.yi, data, enclos = sys.frame(sys.parent()))
    weights <- eval(mf.weights, data, enclos = sys.frame(sys.parent()))
    slab <- eval(mf.slab, data, enclos = sys.frame(sys.parent()))
    subset <- eval(mf.subset, data, enclos = sys.frame(sys.parent()))
    mods <- eval(mf.mods, data, enclos = sys.frame(sys.parent()))
    scale <- eval(mf.scale, data, enclos = sys.frame(sys.parent()))
    ai <- bi <- ci <- di <- x1i <- x2i <- t1i <- t2i <- NA
    is.formula <- FALSE
    if (!is.null(yi)) {
        if (class(yi) == "formula") {
            options(na.action = "na.pass")
            mods <- model.matrix(yi, data = data)
            attr(mods, "assign") <- NULL
            yi <- model.response(model.frame(yi, data = data))
            options(na.action = na.act)
            names(yi) <- NULL
            intercept <- FALSE
            is.formula <- TRUE
        }
        if (is.matrix(yi)) 
            yi <- as.vector(yi)
        k <- length(yi)
        if (measure == "GEN" && !is.null(attr(yi, "measure"))) 
            measure <- attr(yi, "measure")
        attr(yi, "measure") <- measure
        mf.vi <- mf[[match("vi", names(mf))]]
        mf.sei <- mf[[match("sei", names(mf))]]
        mf.ni <- mf[[match("ni", names(mf))]]
        vi <- eval(mf.vi, data, enclos = sys.frame(sys.parent()))
        sei <- eval(mf.sei, data, enclos = sys.frame(sys.parent()))
        ni <- eval(mf.ni, data, enclos = sys.frame(sys.parent()))
        if (is.null(vi)) {
            if (is.null(sei)) {
                stop("Need to specify vi or sei argument.")
            }
            else {
                vi <- sei^2
            }
        }
        if (is.matrix(vi)) 
            vi <- as.vector(vi)
        if (length(vi) == 1) 
            vi <- rep(vi, k)
        if (length(vi) != k) 
            stop("Length of yi and vi (or sei) vectors are not the same.")
        if (is.null(ni) && !is.null(attr(yi, "ni"))) 
            ni <- attr(yi, "ni")
        if (!is.null(ni) && length(ni) != k) 
            ni <- NULL
        if (!is.null(ni)) 
            attr(yi, "ni") <- ni
        if (is.null(slab)) {
            if (!is.null(attr(yi, "slab"))) 
                slab <- attr(yi, "slab")
            if (is.null(slab) && length(slab) != k) 
                slab <- NULL
        }
        if (!is.null(subset)) {
            yi <- yi[subset]
            vi <- vi[subset]
            ni <- ni[subset]
            attr(yi, "measure") <- measure
            attr(yi, "ni") <- ni
        }
    }
    else {
        if (is.element(measure, c("RR", "OR", "PETO", "RD", "AS", 
            "PHI", "YUQ", "YUY", "RTET", "PBIT", "OR2D", "OR2DN", 
            "OR2DL"))) {
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
            dat <- escalc(measure = measure, ai = ai, bi = bi, 
                ci = ci, di = di, add = add, to = to, drop00 = drop00, 
                vtype = vtype)
        }
        if (is.element(measure, c("IRR", "IRD", "IRSD"))) {
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
        if (is.element(measure, c("MD", "SMD", "SMDH", "ROM", 
            "RPB", "RBIS", "D2OR", "D2ORN", "D2ORL"))) {
            mf.m1i <- mf[[match("m1i", names(mf))]]
            mf.m2i <- mf[[match("m2i", names(mf))]]
            mf.sd1i <- mf[[match("sd1i", names(mf))]]
            mf.sd2i <- mf[[match("sd2i", names(mf))]]
            mf.n1i <- mf[[match("n1i", names(mf))]]
            mf.n2i <- mf[[match("n2i", names(mf))]]
            m1i <- eval(mf.m1i, data, enclos = sys.frame(sys.parent()))
            m2i <- eval(mf.m2i, data, enclos = sys.frame(sys.parent()))
            sd1i <- eval(mf.sd1i, data, enclos = sys.frame(sys.parent()))
            sd2i <- eval(mf.sd2i, data, enclos = sys.frame(sys.parent()))
            n1i <- eval(mf.n1i, data, enclos = sys.frame(sys.parent()))
            n2i <- eval(mf.n2i, data, enclos = sys.frame(sys.parent()))
            k <- length(m1i)
            if (!is.null(subset)) {
                m1i <- m1i[subset]
                m2i <- m2i[subset]
                sd1i <- sd1i[subset]
                sd2i <- sd2i[subset]
                n1i <- n1i[subset]
                n2i <- n2i[subset]
            }
            dat <- escalc(measure = measure, m1i = m1i, m2i = m2i, 
                sd1i = sd1i, sd2i = sd2i, n1i = n1i, n2i = n2i, 
                vtype = vtype)
        }
        if (is.element(measure, c("COR", "UCOR", "ZCOR"))) {
            mf.ri <- mf[[match("ri", names(mf))]]
            mf.ni <- mf[[match("ni", names(mf))]]
            ri <- eval(mf.ri, data, enclos = sys.frame(sys.parent()))
            ni <- eval(mf.ni, data, enclos = sys.frame(sys.parent()))
            k <- length(ri)
            if (!is.null(subset)) {
                ri <- ri[subset]
                ni <- ni[subset]
            }
            dat <- escalc(measure = measure, ri = ri, ni = ni, 
                vtype = vtype)
        }
        if (is.element(measure, c("PR", "PLN", "PLO", "PAS", 
            "PFT"))) {
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
            dat <- escalc(measure = measure, xi = xi, mi = mi, 
                add = add, to = to, vtype = vtype)
        }
        if (is.element(measure, c("IR", "IRLN", "IRS", "IRFT"))) {
            mf.xi <- mf[[match("xi", names(mf))]]
            mf.ti <- mf[[match("ti", names(mf))]]
            xi <- eval(mf.xi, data, enclos = sys.frame(sys.parent()))
            ti <- eval(mf.ti, data, enclos = sys.frame(sys.parent()))
            k <- length(xi)
            if (!is.null(subset)) {
                xi <- xi[subset]
                ti <- ti[subset]
            }
            dat <- escalc(measure = measure, xi = xi, ti = ti, 
                add = add, to = to, vtype = vtype)
        }
        if (is.element(measure, c("MN"))) {
            mf.mi <- mf[[match("mi", names(mf))]]
            mf.sdi <- mf[[match("sdi", names(mf))]]
            mf.ni <- mf[[match("ni", names(mf))]]
            mi <- eval(mf.mi, data, enclos = sys.frame(sys.parent()))
            sdi <- eval(mf.sdi, data, enclos = sys.frame(sys.parent()))
            ni <- eval(mf.ni, data, enclos = sys.frame(sys.parent()))
            k <- length(mi)
            if (!is.null(subset)) {
                mi <- mi[subset]
                sdi <- sdi[subset]
                ni <- ni[subset]
            }
            dat <- escalc(measure = measure, mi = mi, sdi = sdi, 
                ni = ni, vtype = vtype)
        }
        if (is.element(measure, c("MC", "SMCC", "SMCR", "SMCRH", 
            "ROMC"))) {
            mf.m1i <- mf[[match("m1i", names(mf))]]
            mf.m2i <- mf[[match("m2i", names(mf))]]
            mf.sd1i <- mf[[match("sd1i", names(mf))]]
            mf.sd2i <- mf[[match("sd2i", names(mf))]]
            mf.ri <- mf[[match("ri", names(mf))]]
            mf.ni <- mf[[match("ni", names(mf))]]
            m1i <- eval(mf.m1i, data, enclos = sys.frame(sys.parent()))
            m2i <- eval(mf.m2i, data, enclos = sys.frame(sys.parent()))
            sd1i <- eval(mf.sd1i, data, enclos = sys.frame(sys.parent()))
            sd2i <- eval(mf.sd2i, data, enclos = sys.frame(sys.parent()))
            ri <- eval(mf.ri, data, enclos = sys.frame(sys.parent()))
            ni <- eval(mf.ni, data, enclos = sys.frame(sys.parent()))
            k <- length(m1i)
            if (!is.null(subset)) {
                m1i <- m1i[subset]
                m2i <- m2i[subset]
                sd1i <- sd1i[subset]
                sd2i <- sd2i[subset]
                ni <- ni[subset]
                ri <- ri[subset]
            }
            dat <- escalc(measure = measure, m1i = m1i, m2i = m2i, 
                sd1i = sd1i, sd2i = sd2i, ri = ri, ni = ni, vtype = vtype)
        }
        if (is.element(measure, c("ARAW", "AHW", "ABT"))) {
            mf.ai <- mf[[match("ai", names(mf))]]
            mf.mi <- mf[[match("mi", names(mf))]]
            mf.ni <- mf[[match("ni", names(mf))]]
            ai <- eval(mf.ai, data, enclos = sys.frame(sys.parent()))
            mi <- eval(mf.mi, data, enclos = sys.frame(sys.parent()))
            ni <- eval(mf.ni, data, enclos = sys.frame(sys.parent()))
            k <- length(ai)
            if (!is.null(subset)) {
                ai <- ai[subset]
                mi <- mi[subset]
                ni <- ni[subset]
            }
            dat <- escalc(measure = measure, ai = ai, mi = mi, 
                ni = ni, vtype = vtype)
        }
        if (is.element(measure, "GEN")) 
            stop("Specify the desired outcome measure via the 'measure' argument.")
        yi <- dat$yi
        vi <- dat$vi
        ni <- attr(yi, "ni")
    }
    if (length(weights) == 1) 
        weights <- rep(weights, k)
    if (!is.null(weights) && (length(weights) != k)) 
        stop("Length of yi and weights vectors are not the same.")
    if (!is.null(subset)) 
        weights <- weights[subset]
    if (very.verbose) 
        message("Creating model matrix ...")
    if (class(mods) == "formula") {
        options(na.action = "na.pass")
        mods <- model.matrix(mods, data = data)
        attr(mods, "assign") <- NULL
        options(na.action = na.act)
        intercept <- FALSE
        is.formula <- TRUE
    }
    if (is.vector(mods)) 
        mods <- cbind(mods)
    if (is.data.frame(mods)) 
        mods <- as.matrix(mods)
    if (is.character(mods)) 
        stop("Model matrix contains character variables.")
    if (!is.null(mods) && (nrow(mods) != k)) 
        stop("Number of rows of the model matrix does not match length of yi argument.")
    if (class(scale) == "formula") {
        options(na.action = "na.pass")
        Z <- model.matrix(scale, data = data)
        attr(Z, "assign") <- NULL
        options(na.action = na.act)
        model <- "rma.tau2"
        if (nrow(Z) != k) 
            stop("Number of rows of the model matrix for tau2 does not match length of yi argument.")
    }
    else {
        Z <- NULL
        model <- "rma.uni"
    }
    ids <- seq_len(k)
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
        Z <- Z[subset, , drop = FALSE]
    }
    if (anyDuplicated(slab)) 
        slab <- make.unique(as.character(slab))
    attr(yi, "slab") <- slab
    k <- length(yi)
    if (any(vi <= 0, na.rm = TRUE)) {
        allvipos <- FALSE
        warning("There are outcomes with non-positive sampling variances.")
        vi.neg <- vi < 0
        if (any(vi.neg, na.rm = TRUE)) {
            vi[vi.neg] <- 0
            warning("Negative sampling variances constrained to zero.")
        }
    }
    else {
        allvipos <- TRUE
    }
    if (any(weights < 0, na.rm = TRUE)) 
        stop("Negative weights not allowed.")
    if (any(is.infinite(weights))) 
        stop("Infinite weights not allowed.")
    ai.f <- ai
    bi.f <- bi
    ci.f <- ci
    di.f <- di
    x1i.f <- x1i
    x2i.f <- x2i
    t1i.f <- t1i
    t2i.f <- t2i
    yi.f <- yi
    vi.f <- vi
    weights.f <- weights
    ni.f <- ni
    mods.f <- mods
    k.f <- k
    YVXZW.na <- is.na(yi) | is.na(vi) | if (is.null(mods)) 
        FALSE
    else apply(is.na(mods), 1, any) | if (is.null(Z)) 
        FALSE
    else apply(is.na(Z), 1, any) | if (is.null(weights)) 
        FALSE
    else is.na(weights)
    if (any(YVXZW.na)) {
        if (very.verbose) 
            message("Handling NAs ...")
        not.na <- !YVXZW.na
        if (na.act == "na.omit" || na.act == "na.exclude" || 
            na.act == "na.pass") {
            yi <- yi[not.na]
            vi <- vi[not.na]
            weights <- weights[not.na]
            ni <- ni[not.na]
            mods <- mods[not.na, , drop = FALSE]
            Z <- Z[not.na, , drop = FALSE]
            k <- length(yi)
            warning("Studies with NAs omitted from model fitting.")
            attr(yi, "measure") <- measure
            attr(yi, "ni") <- ni
        }
        if (na.act == "na.fail") 
            stop("Missing values in data.")
    }
    else {
        not.na <- rep(TRUE, k)
    }
    if (k < 1) 
        stop("Processing terminated since k = 0.")
    if (k == 1) {
        method <- "FE"
        knha <- FALSE
    }
    if (is.null(mods) && !intercept) {
        warning("Must either include an intercept and/or moderators in model.\n  Coerced intercept into the model.")
        intercept <- TRUE
    }
    if (intercept) {
        X <- cbind(intrcpt = rep(1, k), mods)
        X.f <- cbind(intrcpt = rep(1, k.f), mods.f)
    }
    else {
        X <- mods
        X.f <- mods.f
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
    tmp <- lm(yi ~ X - 1)
    coef.na <- is.na(coef(tmp))
    if (any(coef.na)) {
        warning("Redundant predictors dropped from the model.")
        X <- X[, !coef.na, drop = FALSE]
        X.f <- X.f[, !coef.na, drop = FALSE]
    }
    p <- NCOL(X)
    if ((p == 1L) && (all(sapply(X, identical, 1)))) {
        int.only <- TRUE
    }
    else {
        int.only <- FALSE
    }
    if (method == "FE") {
        if (p > k) 
            stop("Number of parameters to be estimated is larger than the number of observations.")
    }
    else {
        if (is.numeric(tau2)) {
            if (p > k) 
                stop("Number of parameters to be estimated is larger than the number of observations.")
        }
        else {
            if ((p + 1) > k) 
                stop("Number of parameters to be estimated is larger than the number of observations.")
        }
    }
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
    con <- list(verbose = FALSE, tau2.init = NULL, tau2.min = 0, 
        tau2.max = 100, threshold = 10^-5, maxiter = 100, stepadj = 1, 
        REMLf = TRUE, tol = 1e-07, optimizer = "nlminb", optmethod = "BFGS")
    con.pos <- pmatch(names(control), names(con))
    con[c(na.omit(con.pos))] <- control[!is.na(con.pos)]
    if (verbose) 
        con$verbose <- verbose
    verbose <- con$verbose
    if (any(eigen(crossprod(X), symmetric = TRUE, only.values = TRUE)$values <= 
        con$tol)) 
        stop("Model matrix not of full rank. Cannot fit model.")
    iter <- 0
    se.tau2 <- I2 <- H2 <- QE <- QEp <- NA
    s2w <- 1
    alpha <- ifelse(level > 1, (100 - level)/100, 1 - level)
    Y <- as.matrix(yi)
    if (model == "rma.tau2") {
        if (method != "ML" && method != "REML") 
            stop("Location-scale model can only be fitted with ML or REML estimation.")
        if (!weighted) 
            stop("Cannot use weighted estimation to fit location-scale model.")
        if (!is.null(weights)) 
            stop("Cannot use user-defined weights for location-scale model.")
        if (any(eigen(crossprod(Z), symmetric = TRUE, only.values = TRUE)$values <= 
            con$tol)) 
            stop("Model matrix for scale part of the model not of full rank. Cannot fit model.")
        optimizer <- match.arg(con$optimizer, c("optim", "nlminb", 
            "uobyqa", "newuoa", "bobyqa", "nloptr"))
        optmethod <- con$optmethod
        optcontrol <- control[is.na(con.pos)]
        if (length(optcontrol) == 0) 
            optcontrol <- list()
        if (optimizer == "nloptr" && !is.element("algorithm", 
            names(optcontrol))) 
            optcontrol$algorithm <- "NLOPT_LN_BOBYQA"
        if (optimizer == "nloptr" && !is.element("ftol_rel", 
            names(optcontrol))) 
            optcontrol$ftol_rel <- 1e-08
        reml <- ifelse(method == "REML", TRUE, FALSE)
        if (is.element(optimizer, c("uobyqa", "newuoa", "bobyqa"))) {
            if (!requireNamespace("minqa", quietly = TRUE)) 
                stop("Please install the 'minqa' package to use this optimizer.")
        }
        if (optimizer == "nloptr") {
            if (!requireNamespace("nloptr", quietly = TRUE)) 
                stop("Please install the 'nloptr' package to use this optimizer.")
        }
        if (!requireNamespace("numDeriv", quietly = TRUE)) 
            stop("Please install the 'numDeriv' package to fit a location-scale model.")
        p.tau2 <- NCOL(Z)
        if (!is.null(con$tau2.init)) {
            if (length(con$tau2.init) != p.tau2) 
                stop(paste("Length of 'tau2.init' argument (", 
                  length(con$tau2.init), ") does not match actual number of parameters (", 
                  p.tau2, ").", sep = ""))
        }
        else {
            con$tau2.init <- rep(1e-04, p.tau2)
        }
        if (optimizer == "optim") 
            par.arg <- "par"
        if (optimizer == "nlminb") 
            par.arg <- "start"
        if (is.element(optimizer, c("uobyqa", "newuoa", "bobyqa"))) {
            par.arg <- "par"
            optimizer <- paste0("minqa::", optimizer)
        }
        if (optimizer == "nloptr") {
            par.arg <- "x0"
            optimizer <- paste0("nloptr::nloptr")
        }
        optcall <- paste(optimizer, "(", par.arg, "=con$tau2.init, .ll.rma.tau2, ", 
            ifelse(optimizer == "optim", "method=optmethod, ", 
                ""), "yi=yi, vi=vi, X=X, Z=Z, reml=reml,\n                                                    k=k, p=p,\n                                                    verbose=verbose, digits=digits, REMLf=con$REMLf, ", 
            ifelse(optimizer == "nloptr::nloptr", "opts=optcontrol)", 
                "control=optcontrol)"), sep = "")
        opt.res <- try(eval(parse(text = optcall)), silent = !verbose)
        if (inherits(opt.res, "try-error")) 
            stop("Error during optimization for location-scale model.")
        if (is.element(optimizer, c("optim", "nlminb")) && opt.res$convergence != 
            0) 
            stop(paste0("Optimizer (", optimizer, ") did not achieve convergence (convergence = ", 
                opt.res$convergence, ")."))
        if (is.element(optimizer, c("minqa::uobyqa", "minqa::newuoa", 
            "minqa::bobyqa")) && opt.res$ierr != 0) 
            stop(paste0("Optimizer (", optimizer, ") did not achieve convergence (ierr = ", 
                opt.res$ierr, ")."))
        if (optimizer == "nloptr::nloptr" && !(opt.res$status >= 
            1 && opt.res$status <= 4)) 
            stop(paste0("Optimizer (", optimizer, ") did not achieve convergence (status = ", 
                opt.res$status, ")."))
        if (verbose) {
            cat("\n")
            print(opt.res)
        }
        if (optimizer == "nloptr::nloptr") 
            opt.res$par <- opt.res$solution
        vb.tau2 <- matrix(NA, nrow = p.tau2, ncol = p.tau2)
        se.tau2 <- rep(NA, p.tau2)
        h <- try(numDeriv::hessian(.ll.rma.tau2, x = opt.res$par, 
            yi = yi, vi = vi, X = X, Z = Z, reml = reml, k = k, 
            p = p, verbose = FALSE, digits = digits, REMLf = con$REMLf), 
            silent = TRUE)
        if (!inherits(h, "try-error")) {
            chol.h <- try(chol(h), silent = TRUE)
            if (!inherits(chol.h, "try-error")) {
                vb.tau2 <- chol2inv(chol.h)
                se.tau2 <- sqrt(diag(vb.tau2))
            }
        }
        b.tau2 <- cbind(opt.res$par)
        colnames(Z)[grep("(Intercept)", colnames(Z))] <- "intrcpt"
        rownames(b.tau2) <- rownames(vb.tau2) <- colnames(vb.tau2) <- colnames(Z)
        crit <- qnorm(alpha/2, lower.tail = FALSE)
        names(se.tau2) <- NULL
        zval.tau2 <- c(b.tau2/se.tau2)
        pval.tau2 <- 2 * pnorm(abs(zval.tau2), lower.tail = FALSE)
        ci.lb.tau2 <- c(b.tau2 - crit * se.tau2)
        ci.ub.tau2 <- c(b.tau2 + crit * se.tau2)
        tau2 <- exp(as.vector(Z %*% b.tau2))
        tau2.fix <- FALSE
    }
    else {
        if (is.numeric(tau2)) {
            tau2.fix <- TRUE
            tau2.val <- tau2
        }
        else {
            tau2.fix <- FALSE
            tau2.val <- NA
        }
        if (very.verbose && !tau2.fix) 
            message("Estimating tau^2 value ...")
        if (method == "HS") {
            if (!allvipos) 
                stop("HS estimator cannot be used with non-positive sampling variances.")
            wi <- 1/vi
            W <- diag(wi, nrow = k, ncol = k)
            stXWX <- .invcalc(X = X, W = W, k = k)
            P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
            RSS <- crossprod(Y, P) %*% Y
            tau2 <- ifelse(tau2.fix, tau2.val, (RSS - k)/sum(wi))
        }
        if (is.element(method, c("HE", "ML", "REML", "EB"))) {
            stXX <- .invcalc(X = X, W = diag(k), k = k)
            P <- diag(k) - X %*% tcrossprod(stXX, X)
            RSS <- crossprod(Y, P) %*% Y
            V <- diag(vi, nrow = k, ncol = k)
            PV <- P %*% V
            trPV <- .tr(PV)
            tau2 <- ifelse(tau2.fix, tau2.val, (RSS - trPV)/(k - 
                p))
        }
        if (method == "DL") {
            if (!allvipos) 
                stop("DL estimator cannot be used with non-positive sampling variances.")
            wi <- 1/vi
            W <- diag(wi, nrow = k, ncol = k)
            stXWX <- .invcalc(X = X, W = W, k = k)
            P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
            RSS <- crossprod(Y, P) %*% Y
            trP <- .tr(P)
            tau2 <- ifelse(tau2.fix, tau2.val, (RSS - (k - p))/trP)
        }
        if (method == "GENQ") {
            if (is.null(weights)) 
                stop("Must specify 'weights' when method='GENQ'.")
            A <- diag(weights, nrow = k, ncol = k)
            stXAX <- .invcalc(X = X, W = A, k = k)
            P <- A - A %*% X %*% stXAX %*% t(X) %*% A
            RSS <- crossprod(Y, P) %*% Y
            V <- diag(vi, nrow = k, ncol = k)
            PV <- P %*% V
            trP <- .tr(P)
            trPV <- .tr(PV)
            tau2 <- ifelse(tau2.fix, tau2.val, (RSS - trPV)/trP)
        }
        if (method == "SJ") {
            if (is.null(con$tau2.init)) {
                tau2.0 <- c(var(yi) * (k - 1)/k)
            }
            else {
                tau2.0 <- con$tau2.init
            }
            wi <- 1/(vi + tau2.0)
            W <- diag(wi, nrow = k, ncol = k)
            stXWX <- .invcalc(X = X, W = W, k = k)
            P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
            RSS <- crossprod(Y, P) %*% Y
            V <- diag(vi, nrow = k, ncol = k)
            PV <- P %*% V
            tau2 <- ifelse(tau2.fix, tau2.val, tau2.0 * RSS/(k - 
                p))
        }
        if (method == "DLIT") {
            conv <- 1
            change <- con$threshold + 1
            if (is.null(con$tau2.init)) {
                tau2 <- 0
            }
            else {
                tau2 <- con$tau2.init
            }
            while (change > con$threshold) {
                if (verbose) 
                  cat("Iteration", iter, "\ttau^2 =", formatC(tau2, 
                    format = "f", digits = digits), "\n")
                iter <- iter + 1
                tau2.old <- tau2
                wi <- 1/(vi + tau2)
                if (any(is.infinite(wi))) 
                  stop("Division by zero when computing the inverse variance weights.")
                W <- diag(wi, nrow = k, ncol = k)
                stXWX <- .invcalc(X = X, W = W, k = k)
                P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
                RSS <- crossprod(Y, P) %*% Y
                trP <- .tr(P)
                tau2 <- ifelse(tau2.fix, tau2.val, (RSS - (k - 
                  p))/trP)
                tau2[tau2 < con$tau2.min] <- con$tau2.min
                change <- abs(tau2.old - tau2)
                if (iter > con$maxiter) {
                  conv <- 0
                  break
                }
            }
            if (conv == 0L) 
                stop("Algorithm did not converge.")
        }
        if (method == "SJIT") {
            conv <- 1
            change <- con$threshold + 1
            if (is.null(con$tau2.init)) {
                tau2 <- var(yi) * (k - 1)/k
                tau2.0 <- tau2
            }
            else {
                tau2 <- con$tau2.init
                tau2.0 <- tau2
            }
            while (change > con$threshold) {
                if (verbose) 
                  cat("Iteration", iter, "\ttau^2 =", formatC(tau2, 
                    format = "f", digits = digits), "\n")
                iter <- iter + 1
                tau2.old <- tau2
                wi <- 1/(vi + tau2)
                W <- diag(wi, nrow = k, ncol = k)
                stXWX <- .invcalc(X = X, W = W, k = k)
                P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
                RSS <- crossprod(Y, P) %*% Y
                V <- diag(vi, nrow = k, ncol = k)
                PV <- P %*% V
                tau2 <- ifelse(tau2.fix, tau2.val, tau2 * RSS/(k - 
                  p))
                change <- abs(tau2.old - tau2)
                if (iter > con$maxiter) {
                  conv <- 0
                  break
                }
            }
            if (conv == 0L) 
                stop("Algorithm did not converge.")
        }
        if (method == "PM") {
            if (!allvipos) 
                stop("PM estimator cannot be used with non-positive sampling variances.")
            if (.QE.func(con$tau2.min, Y = Y, vi = vi, X = X, 
                k = k, objective = k - p) < con$tau2.min) {
                tau2 <- con$tau2.min
            }
            else {
                tau2 <- ifelse(tau2.fix, tau2.val, try(uniroot(.QE.func, 
                  interval = c(con$tau2.min, con$tau2.max), tol = con$threshold, 
                  maxiter = con$maxiter, Y = Y, vi = vi, X = X, 
                  k = k, objective = k - p, verbose = verbose, 
                  digits = digits, extendInt = "upX")$root, silent = TRUE))
                if (!is.numeric(tau2)) 
                  stop("Error in iterative search for tau2. Try increasing tau2.max or switch to another 'method'.")
            }
            wi <- 1/(vi + tau2)
            W <- diag(wi, nrow = k, ncol = k)
            stXWX <- .invcalc(X = X, W = W, k = k)
            P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
        }
        if (is.element(method, c("ML", "REML", "EB"))) {
            conv <- 1
            change <- con$threshold + 1
            if (is.null(con$tau2.init)) {
                tau2 <- max(0, tau2, con$tau2.min)
            }
            else {
                tau2 <- con$tau2.init
            }
            while (change > con$threshold) {
                if (verbose) 
                  cat("Iteration", iter, "\ttau^2 =", formatC(tau2, 
                    format = "f", digits = digits), "\n")
                iter <- iter + 1
                tau2.old <- tau2
                wi <- 1/(vi + tau2)
                if (any(is.infinite(wi))) 
                  stop("Division by zero when computing the inverse variance weights.")
                W <- diag(wi, nrow = k, ncol = k)
                stXWX <- .invcalc(X = X, W = W, k = k)
                P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
                if (method == "ML") {
                  PP <- P %*% P
                  adj <- (crossprod(Y, PP) %*% Y - sum(wi))/sum(wi^2)
                }
                if (method == "REML") {
                  PP <- P %*% P
                  adj <- (crossprod(Y, PP) %*% Y - .tr(P))/.tr(PP)
                }
                if (method == "EB") {
                  adj <- (crossprod(Y, P) %*% Y * k/(k - p) - 
                    k)/sum(wi)
                }
                adj <- adj * con$stepadj
                while (tau2 + adj < con$tau2.min) {
                  adj <- adj/2
                }
                tau2 <- ifelse(tau2.fix, tau2.val, tau2 + adj)
                change <- abs(tau2.old - tau2)
                if (iter > con$maxiter) {
                  conv <- 0
                  break
                }
            }
            if (conv == 0L) 
                stop("Fisher scoring algorithm did not converge. See 'help(rma)' for possible remedies.")
        }
        tau2 <- max(con$tau2.min, c(tau2))
        if (verbose && is.element(method, c("ML", "REML", "EB"))) {
            cat("Iteration", iter, "\ttau^2 =", formatC(tau2, 
                format = "f", digits = digits), "\n")
            cat("Fisher scoring algorithm converged after", iter, 
                "iterations.\n")
        }
        if (method == "HS") {
            se.tau2 <- sqrt(1/sum(wi)^2 * (2 * (k - p) + 4 * 
                max(tau2, 0) * .tr(P) + 2 * max(tau2, 0)^2 * 
                sum(P * P)))
        }
        if (method == "HE") {
            se.tau2 <- sqrt(1/(k - p)^2 * (2 * sum(PV * t(PV)) + 
                4 * max(tau2, 0) * trPV + 2 * max(tau2, 0)^2 * 
                (k - p)))
        }
        if (method == "DL" || method == "DLIT") {
            se.tau2 <- sqrt(1/trP^2 * (2 * (k - p) + 4 * max(tau2, 
                0) * trP + 2 * max(tau2, 0)^2 * sum(P * P)))
        }
        if (method == "GENQ") {
            se.tau2 <- sqrt(1/trP^2 * (2 * sum(PV * t(PV)) + 
                4 * max(tau2, 0) * sum(PV * P) + 2 * max(tau2, 
                0)^2 * sum(P * P)))
        }
        if (method == "SJ") {
            se.tau2 <- sqrt(tau2.0^2/(k - p)^2 * (2 * sum(PV * 
                t(PV)) + 4 * max(tau2, 0) * sum(PV * P) + 2 * 
                max(tau2, 0)^2 * sum(P * P)))
        }
        if (method == "ML") {
            se.tau2 <- sqrt(2/sum(wi^2))
        }
        if (method == "REML") {
            se.tau2 <- sqrt(2/sum(P * P))
        }
        if (method == "EB" || method == "PM" || method == "SJIT") {
            V <- diag(vi, nrow = k, ncol = k)
            PV <- P %*% V
            se.tau2 <- sqrt(2 * k^2/(k - p)/sum(wi)^2)
        }
    }
    if (method == "FE") 
        tau2 <- 0
    if (very.verbose) 
        message("Model fitting ...")
    wi <- 1/(vi + tau2)
    W <- diag(wi, nrow = k, ncol = k)
    M <- diag(vi + tau2, nrow = k, ncol = k)
    if (weighted) {
        if (is.null(weights)) {
            if (any(is.infinite(wi))) 
                stop("Division by zero when computing the inverse variance weights.")
            stXWX <- .invcalc(X = X, W = W, k = k)
            b <- stXWX %*% crossprod(X, W) %*% Y
            vb <- stXWX
            RSS.f <- sum(wi * (yi - X %*% b)^2)
        }
        else {
            A <- diag(weights, nrow = k, ncol = k)
            stXAX <- .invcalc(X = X, W = A, k = k)
            b <- stXAX %*% crossprod(X, A) %*% Y
            vb <- stXAX %*% t(X) %*% A %*% M %*% A %*% X %*% 
                stXAX
            RSS.f <- sum(wi * (yi - X %*% b)^2)
        }
        if (knha) {
            if (RSS.f <= .Machine$double.eps) {
                s2w <- 1
            }
            else {
                s2w <- RSS.f/(k - p)
            }
            vb <- s2w * vb
            if (method == "FE") 
                warning("Knapp & Hartung (2003) method is not meant to be used in the context of fixed-effects models.")
        }
        QM <- c(t(b)[btt] %*% chol2inv(chol(vb[btt, btt])) %*% 
            b[btt])
    }
    else {
        stXX <- .invcalc(X = X, W = diag(k), k = k)
        b <- stXX %*% crossprod(X, Y)
        vb <- tcrossprod(stXX, X) %*% M %*% X %*% stXX
        RSS.f <- sum(wi * (yi - X %*% b)^2)
        if (knha) {
            stXWX <- .invcalc(X = X, W = W, k = k)
            b.knha <- stXWX %*% crossprod(X, W) %*% Y
            RSS.knha <- sum(wi * (yi - X %*% b.knha)^2)
            if (RSS.knha <= .Machine$double.eps) {
                s2w <- 1
            }
            else {
                s2w <- RSS.knha/(k - p)
            }
            vb <- s2w * vb
            if (method == "FE") 
                warning("Knapp & Hartung (2003) method is not meant to be used in the context of fixed-effects models.")
        }
        QM <- c(t(b)[btt] %*% chol2inv(chol(vb[btt, btt])) %*% 
            b[btt])
    }
    rownames(b) <- rownames(vb) <- colnames(vb) <- colnames(X)
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
        message("Heterogeneity testing ...")
    if (allvipos) {
        wi <- 1/vi
        W.FE <- diag(wi, nrow = k, ncol = k)
        stXWX <- .invcalc(X = X, W = W.FE, k = k)
        P <- W.FE - W.FE %*% X %*% stXWX %*% crossprod(X, W.FE)
        QE <- max(0, c(crossprod(Y, P) %*% Y))
        QEp <- ifelse(k - p >= 1, pchisq(QE, df = k - p, lower.tail = FALSE), 
            1)
        vi.avg <- (k - p)/.tr(P)
        I2 <- 100 * mean(tau2)/(vi.avg + mean(tau2))
        H2 <- mean(tau2)/vi.avg + 1
    }
    else {
        warning(paste0("Cannot compute ", ifelse(int.only, "Q", 
            "QE"), "-test, I^2, or H^2 with non-positive sampling variances."))
    }
    if (!int.only && int.incl && method != "FE" && model == "rma.uni") {
        if (very.verbose) {
            message("Fitting RE model for R^2 computation ...")
            res.RE <- try(rma(yi, vi, weights = weights, method = method, 
                weighted = weighted, knha = knha, verbose = ifelse(verbose, 
                  TRUE, FALSE), control = con, digits = digits), 
                silent = FALSE)
        }
        else {
            res.RE <- suppressWarnings(try(rma(yi, vi, weights = weights, 
                method = method, weighted = weighted, knha = knha, 
                verbose = ifelse(verbose, TRUE, FALSE), control = con, 
                digits = digits), silent = FALSE))
        }
        if (!inherits(res.RE, "try-error")) {
            tau2.RE <- res.RE$tau2
            if (identical(tau2.RE, 0)) {
                R2 <- NA
            }
            else {
                R2 <- round(max(0, 100 * (tau2.RE - tau2)/tau2.RE), 
                  2)
            }
        }
        else {
            R2 <- NA
        }
    }
    else {
        R2 <- NULL
    }
    if (very.verbose) 
        message("Computing fit statistics and log likelihood ...")
    parms <- p + ifelse(model == "rma.uni", ifelse(method == 
        "FE" || tau2.fix, 0, 1), p.tau2)
    ll.ML <- -1/2 * (k) * log(2 * base::pi) - 1/2 * sum(log(vi + 
        tau2)) - 1/2 * RSS.f
    ll.REML <- -1/2 * (k - p) * log(2 * base::pi) + ifelse(con$REMLf, 
        1/2 * determinant(crossprod(X), logarithm = TRUE)$modulus, 
        0) - 1/2 * sum(log(vi + tau2)) - 1/2 * determinant(crossprod(X, 
        W) %*% X, logarithm = TRUE)$modulus - 1/2 * RSS.f
    if (k - p > 0L) {
        dev.ML <- -2 * (ll.ML - sum(dnorm(yi, mean = yi, sd = sqrt(vi), 
            log = TRUE)))
    }
    else {
        dev.ML <- 0
    }
    AIC.ML <- -2 * ll.ML + 2 * parms
    BIC.ML <- -2 * ll.ML + parms * log(k)
    AICc.ML <- -2 * ll.ML + 2 * parms * max(k, parms + 2)/(max(k, 
        parms + 2) - parms - 1)
    dev.REML <- -2 * (ll.REML - 0)
    AIC.REML <- -2 * ll.REML + 2 * parms
    BIC.REML <- -2 * ll.REML + parms * log(k - p)
    AICc.REML <- -2 * ll.REML + 2 * parms * max(k - p, parms + 
        2)/(max(k - p, parms + 2) - parms - 1)
    fit.stats <- matrix(c(ll.ML, dev.ML, AIC.ML, BIC.ML, AICc.ML, 
        ll.REML, dev.REML, AIC.REML, BIC.REML, AICc.REML), ncol = 2, 
        byrow = FALSE)
    dimnames(fit.stats) <- list(c("ll", "dev", "AIC", "BIC", 
        "AICc"), c("ML", "REML"))
    fit.stats <- data.frame(fit.stats)
    if (very.verbose) 
        message("Preparing output ...")
    p.eff <- p
    k.eff <- k
    res <- list(b = b, se = se, zval = zval, pval = pval, ci.lb = ci.lb, 
        ci.ub = ci.ub, vb = vb, tau2 = tau2, se.tau2 = se.tau2, 
        tau2.fix = tau2.fix, k = k, k.f = k.f, k.eff = k.eff, 
        p = p, p.eff = p.eff, parms = parms, m = m, QE = QE, 
        QEp = QEp, QM = QM, QMp = QMp, I2 = I2, H2 = H2, R2 = R2, 
        int.only = int.only, int.incl = int.incl, allvipos = allvipos, 
        coef.na = coef.na, yi = yi, vi = vi, X = X, weights = weights, 
        yi.f = yi.f, vi.f = vi.f, X.f = X.f, weights.f = weights.f, 
        M = M, ai.f = ai.f, bi.f = bi.f, ci.f = ci.f, di.f = di.f, 
        x1i.f = x1i.f, x2i.f = x2i.f, t1i.f = t1i.f, t2i.f = t2i.f, 
        ni = ni, ni.f = ni.f, ids = ids, not.na = not.na, subset = subset, 
        slab = slab, slab.null = slab.null, measure = measure, 
        method = method, weighted = weighted, knha = knha, dfs = dfs, 
        s2w = s2w, btt = btt, intercept = intercept, digits = digits, 
        level = level, sparse = FALSE, control = control, verbose = verbose, 
        add = add, to = to, drop00 = drop00, fit.stats = fit.stats, 
        version = packageVersion("metafor"), model = model, call = mf)
    if (model == "rma.tau2") {
        res$b.tau2 <- b.tau2
        res$vb.tau2 <- vb.tau2
        res$se.tau2 <- se.tau2
        res$zval.tau2 <- zval.tau2
        res$pval.tau2 <- pval.tau2
        res$ci.lb.tau2 <- ci.lb.tau2
        res$ci.ub.tau2 <- ci.ub.tau2
        res$Z <- Z
    }
    class(res) <- c("rma.uni", "rma")
    return(res)
}
