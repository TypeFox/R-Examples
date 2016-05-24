escalc.default <-
function (measure, formula, ai, bi, ci, di, n1i, n2i, x1i, x2i, 
    t1i, t2i, m1i, m2i, sd1i, sd2i, xi, mi, ri, ti, sdi, ni, 
    data, slab, subset, add = 1/2, to = "only0", drop00 = FALSE, 
    vtype = "LS", var.names = c("yi", "vi"), add.measure = FALSE, 
    append = TRUE, replace = TRUE, digits = 4, ...) 
{
    if (!is.element(measure, c("RR", "OR", "PETO", "RD", "AS", 
        "PHI", "YUQ", "YUY", "RTET", "PBIT", "OR2D", "OR2DN", 
        "OR2DL", "IRR", "IRD", "IRSD", "MD", "SMD", "SMDH", "ROM", 
        "RPB", "RBIS", "D2OR", "D2ORN", "D2ORL", "COR", "UCOR", 
        "ZCOR", "PR", "PLN", "PLO", "PAS", "PFT", "IR", "IRLN", 
        "IRS", "IRFT", "MN", "MC", "SMCC", "SMCR", "SMCRH", "ROMC", 
        "ARAW", "AHW", "ABT"))) 
        stop("Unknown 'measure' specified.")
    if (!is.element(to, c("all", "only0", "if0all", "none"))) 
        stop("Unknown 'to' argument specified.")
    if (any(!is.element(vtype, c("UB", "LS", "HO", "ST", "CS")), 
        na.rm = TRUE)) 
        stop("Unknown 'vtype' argument specified.")
    if (add.measure) {
        if (length(var.names) == 2) 
            var.names <- c(var.names, "measure")
        if (length(var.names) != 3) 
            stop("Argument var.names must be of length 2 or 3.")
        if (any(var.names != make.names(var.names, unique = TRUE))) {
            var.names <- make.names(var.names, unique = TRUE)
            warning(paste0("Argument 'var.names' does not contain syntactically valid variable names.\n  Variable names adjusted to: var.names = c('", 
                var.names[1], "', '", var.names[2], "', '", var.names[3], 
                "')."))
        }
    }
    else {
        if (length(var.names) == 2) 
            var.names <- var.names[1:2]
        if (length(var.names) != 2) 
            stop("Argument var.names must be of length 2.")
        if (any(var.names != make.names(var.names, unique = TRUE))) {
            var.names <- make.names(var.names, unique = TRUE)
            warning(paste0("Argument 'var.names' does not contain syntactically valid variable names.\n  Variable names adjusted to: var.names = c('", 
                var.names[1], "', '", var.names[2], "')."))
        }
    }
    if (missing(data)) 
        data <- NULL
    no.data <- is.null(data)
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
    if (is.element(measure, c("RR", "OR", "RD", "AS", "PETO", 
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
        if (!is.null(subset)) {
            ai <- ai[subset]
            bi <- bi[subset]
            ci <- ci[subset]
            di <- di[subset]
        }
        if (length(ai) == 0L || length(bi) == 0L || length(ci) == 
            0L || length(di) == 0L) 
            stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments.")
        if (!all(length(ai) == c(length(ai), length(bi), length(ci), 
            length(di)))) 
            stop("Supplied data vectors are not all of the same length.")
        if (any(c(ai, bi, ci, di) < 0, na.rm = TRUE)) 
            stop("One or more counts are negative.")
        ni.u <- ai + bi + ci + di
        k <- length(ai)
        if (drop00) {
            id00 <- c(ai == 0L & ci == 0L) | c(bi == 0L & di == 
                0L)
            id00[is.na(id00)] <- FALSE
            ai[id00] <- NA
            bi[id00] <- NA
            ci[id00] <- NA
            di[id00] <- NA
        }
        if (to == "all") {
            ai <- ai + add
            ci <- ci + add
            bi <- bi + add
            di <- di + add
        }
        if (to == "only0") {
            id0 <- c(ai == 0L | ci == 0L | bi == 0L | di == 0L)
            id0[is.na(id0)] <- FALSE
            ai[id0] <- ai[id0] + add
            ci[id0] <- ci[id0] + add
            bi[id0] <- bi[id0] + add
            di[id0] <- di[id0] + add
        }
        if (to == "if0all") {
            id0 <- c(ai == 0L | ci == 0L | bi == 0L | di == 0L)
            id0[is.na(id0)] <- FALSE
            if (any(id0)) {
                ai <- ai + add
                ci <- ci + add
                bi <- bi + add
                di <- di + add
            }
        }
        n1i <- ai + bi
        n2i <- ci + di
        ni <- n1i + n2i
        p1i <- ai/n1i
        p2i <- ci/n2i
        if (measure == "RR") {
            yi <- log(p1i) - log(p2i)
            vi <- 1/ai - 1/n1i + 1/ci - 1/n2i
        }
        if (is.element(measure, c("OR", "OR2D", "OR2DN", "OR2DL"))) {
            yi <- log(p1i/(1 - p1i)) - log(p2i/(1 - p2i))
            vi <- 1/ai + 1/bi + 1/ci + 1/di
        }
        if (measure == "PETO") {
            xt <- ai + ci
            yt <- bi + di
            Oi <- ai
            Ei <- xt * n1i/ni
            Vi <- xt * yt * (n1i/ni) * (n2i/ni)/(ni - 1)
            yi <- (Oi - Ei)/Vi
            vi <- 1/Vi
        }
        if (measure == "RD") {
            yi <- p1i - p2i
            if (length(vtype) == 1L) 
                vtype <- rep(vtype, k)
            vi <- rep(NA_real_, k)
            mnwp1i <- sum(ai, na.rm = TRUE)/sum(n1i, na.rm = TRUE)
            mnwp2i <- sum(ci, na.rm = TRUE)/sum(n2i, na.rm = TRUE)
            for (i in seq_len(k)) {
                if (vtype[i] == "UB") 
                  vi[i] <- p1i[i] * (1 - p1i[i])/(n1i[i] - 1) + 
                    p2i[i] * (1 - p2i[i])/(n2i[i] - 1)
                if (vtype[i] == "LS") 
                  vi[i] <- p1i[i] * (1 - p1i[i])/n1i[i] + p2i[i] * 
                    (1 - p2i[i])/n2i[i]
                if (vtype[i] == "HO") 
                  vi[i] <- mnwp1i * (1 - mnwp1i)/n1i[i] + mnwp2i * 
                    (1 - mnwp2i)/n2i[i]
            }
        }
        if (measure == "AS") {
            yi <- asin(sqrt(p1i)) - asin(sqrt(p2i))
            vi <- 1/(4 * n1i) + 1/(4 * n2i)
        }
        if (measure == "PHI") {
            yi <- (ai * di - bi * ci)/sqrt((ai + bi) * (ci + 
                di) * (ai + ci) * (bi + di))
            p1. <- (ai + bi)/ni
            p2. <- (ci + di)/ni
            p.1 <- (ai + ci)/ni
            p.2 <- (bi + di)/ni
            vi <- 1/ni * (1 - yi^2 + yi * (1 + 1/2 * yi^2) * 
                (p1. - p2.) * (p.1 - p.2)/sqrt(p1. * p2. * p.1 * 
                p.2) - 3/4 * yi^2 * ((p1. - p2.)^2/(p1. * p2.) + 
                (p.1 - p.2)^2/(p.1 * p.2)))
        }
        if (measure == "YUQ") {
            yi <- (ai/bi)/(ci/di)
            yi <- (yi - 1)/(yi + 1)
            vi <- 1/4 * (1 - yi^2)^2 * (1/ai + 1/bi + 1/ci + 
                1/di)
        }
        if (measure == "YUY") {
            yi <- (ai/bi)/(ci/di)
            yi <- (sqrt(yi) - 1)/(sqrt(yi) + 1)
            vi <- 1/16 * (1 - yi^2)^2 * (1/ai + 1/bi + 1/ci + 
                1/di)
        }
        if (measure == "RTET") {
            yi <- rep(NA_real_, k)
            vi <- rep(NA_real_, k)
            for (i in seq_len(k)) {
                if (is.na(ai[i]) || is.na(bi[i]) || is.na(ci[i]) || 
                  is.na(di[i])) 
                  next
                res <- .rtet(ai[i], bi[i], ci[i], di[i], maxcor = 0.9999)
                if (inherits(res, "try-error")) {
                  next
                }
                else {
                  yi[i] <- res$yi
                  vi[i] <- res$vi
                }
            }
        }
        if (measure == "PBIT") {
            z1i <- qnorm(p1i)
            z2i <- qnorm(p2i)
            yi <- z1i - z2i
            vi <- 2 * pi * p1i * (1 - p1i) * exp(z1i^2)/n1i + 
                2 * pi * p2i * (1 - p2i) * exp(z2i^2)/n2i
        }
        if (is.element(measure, c("OR2D", "OR2DL"))) {
            yi <- sqrt(3)/pi * yi
            vi <- 3/pi^2 * vi
        }
        if (measure == "OR2DN") {
            yi <- yi/1.65
            vi <- vi/1.65^2
        }
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
        if (!is.null(subset)) {
            x1i <- x1i[subset]
            x2i <- x2i[subset]
            t1i <- t1i[subset]
            t2i <- t2i[subset]
        }
        if (length(x1i) == 0L || length(x2i) == 0L || length(t1i) == 
            0L || length(t2i) == 0L) 
            stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments.")
        if (!all(length(x1i) == c(length(x1i), length(x2i), length(t1i), 
            length(t2i)))) 
            stop("Supplied data vectors are not all of the same length.")
        if (any(c(x1i, x2i) < 0, na.rm = TRUE)) 
            stop("One or more counts are negative.")
        if (any(c(t1i, t2i) < 0, na.rm = TRUE)) 
            stop("One or more person-times are negative.")
        ni.u <- t1i + t2i
        if (drop00) {
            id00 <- c(x1i == 0L & x2i == 0L)
            id00[is.na(id00)] <- FALSE
            x1i[id00] <- NA
            x2i[id00] <- NA
        }
        if (to == "all") {
            x1i <- x1i + add
            x2i <- x2i + add
        }
        if (to == "only0") {
            id0 <- c(x1i == 0L | x2i == 0L)
            id0[is.na(id0)] <- FALSE
            x1i[id0] <- x1i[id0] + add
            x2i[id0] <- x2i[id0] + add
        }
        if (to == "if0all") {
            id0 <- c(x1i == 0L | x2i == 0L)
            id0[is.na(id0)] <- FALSE
            if (any(id0)) {
                x1i <- x1i + add
                x2i <- x2i + add
            }
        }
        ir1i <- x1i/t1i
        ir2i <- x2i/t2i
        if (measure == "IRR") {
            yi <- log(ir1i) - log(ir2i)
            vi <- 1/x1i + 1/x2i
        }
        if (measure == "IRD") {
            yi <- ir1i - ir2i
            vi <- ir1i/t1i + ir2i/t2i
        }
        if (measure == "IRSD") {
            yi <- sqrt(ir1i) - sqrt(ir2i)
            vi <- 1/(4 * t1i) + 1/(4 * t2i)
        }
    }
    if (is.element(measure, c("MD", "SMD", "SMDH", "ROM", "RPB", 
        "RBIS", "D2OR", "D2ORN", "D2ORL"))) {
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
        if (!is.null(subset)) {
            m1i <- m1i[subset]
            m2i <- m2i[subset]
            sd1i <- sd1i[subset]
            sd2i <- sd2i[subset]
            n1i <- n1i[subset]
            n2i <- n2i[subset]
        }
        if (length(m1i) == 0L || length(m2i) == 0L || length(sd1i) == 
            0L || length(sd2i) == 0L || length(n1i) == 0L || 
            length(n2i) == 0L) 
            stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments.")
        if (!all(length(m1i) == c(length(m1i), length(m2i), length(sd1i), 
            length(sd2i), length(n1i), length(n2i)))) 
            stop("Supplied data vectors are not all of the same length.")
        if (any(c(sd1i, sd2i) < 0, na.rm = TRUE)) 
            stop("One or more standard deviations are negative.")
        if (any(c(n1i, n2i) < 0, na.rm = TRUE)) 
            stop("One or more sample sizes are negative.")
        ni.u <- n1i + n2i
        k <- length(m1i)
        ni <- ni.u
        mi <- ni - 2
        sdpi <- sqrt(((n1i - 1) * sd1i^2 + (n2i - 1) * sd2i^2)/mi)
        di <- (m1i - m2i)/sdpi
        if (measure == "MD") {
            yi <- m1i - m2i
            if (length(vtype) == 1L) 
                vtype <- rep(vtype, k)
            vi <- rep(NA_real_, k)
            for (i in seq_len(k)) {
                if (vtype[i] == "UB" || vtype[i] == "LS") 
                  vi[i] <- sd1i[i]^2/n1i[i] + sd2i[i]^2/n2i[i]
                if (vtype[i] == "HO") 
                  vi[i] <- sdpi[i]^2 * (1/n1i[i] + 1/n2i[i])
            }
        }
        if (measure == "SMD") {
            warn.before <- getOption("warn")
            options(warn = -1)
            cmi <- .cmicalc(mi)
            options(warn = warn.before)
            yi <- cmi * di
            if (length(vtype) == 1L) 
                vtype <- rep(vtype, k)
            vi <- rep(NA_real_, k)
            mnwyi <- sum(ni * yi, na.rm = TRUE)/sum(ni, na.rm = TRUE)
            for (i in seq_len(k)) {
                if (vtype[i] == "UB") 
                  vi[i] <- 1/n1i[i] + 1/n2i[i] + (1 - (mi[i] - 
                    2)/(mi[i] * cmi[i]^2)) * yi[i]^2
                if (vtype[i] == "LS") 
                  vi[i] <- 1/n1i[i] + 1/n2i[i] + yi[i]^2/(2 * 
                    ni[i])
                if (vtype[i] == "HO") 
                  vi[i] <- 1/n1i[i] + 1/n2i[i] + mnwyi^2/(2 * 
                    ni[i])
            }
        }
        if (measure == "SMDH") {
            warn.before <- getOption("warn")
            options(warn = -1)
            cmi <- .cmicalc(mi)
            options(warn = warn.before)
            si <- sqrt((sd1i^2 + sd2i^2)/2)
            yi <- cmi * (m1i - m2i)/si
            vi <- yi^2 * (sd1i^4/(n1i - 1) + sd2i^4/(n2i - 1))/(2 * 
                (sd1i^2 + sd2i^2)^2) + (sd1i^2/(n1i - 1) + sd2i^2/(n2i - 
                1))/((sd1i^2 + sd2i^2)/2)
            vi <- cmi^2 * vi
        }
        if (measure == "ROM") {
            yi <- log(m1i/m2i)
            if (length(vtype) == 1L) 
                vtype <- rep(vtype, k)
            vi <- rep(NA_real_, k)
            for (i in seq_len(k)) {
                if (vtype[i] == "LS") 
                  vi[i] <- sd1i[i]^2/(n1i[i] * m1i[i]^2) + sd2i[i]^2/(n2i[i] * 
                    m2i[i]^2)
                if (vtype[i] == "HO") 
                  vi[i] <- sdpi[i]^2/(n1i[i] * m1i[i]^2) + sdpi[i]^2/(n2i[i] * 
                    m2i[i]^2)
            }
        }
        if (is.element(measure, c("RPB", "RBIS"))) {
            hi <- mi/n1i + mi/n2i
            yi <- di/sqrt(di^2 + hi)
            if (measure == "RPB") {
                if (length(vtype) == 1L) 
                  vtype <- rep(vtype, k)
                vi <- rep(NA_real_, k)
                for (i in seq_len(k)) {
                  if (vtype[i] == "ST" || vtype[i] == "LS") 
                    vi[i] <- hi[i]^2/(hi[i] + di[i]^2)^3 * (1/n1i[i] + 
                      1/n2i[i] + di[i]^2/(2 * ni[i]))
                  if (vtype[i] == "CS") 
                    vi[i] <- (1 - yi[i]^2)^2 * (ni[i] * yi[i]^2/(4 * 
                      n1i[i] * n2i[i]) + (2 - 3 * yi[i]^2)/(2 * 
                      ni[i]))
                }
            }
        }
        if (measure == "RBIS") {
            p1i <- n1i/ni
            p2i <- n2i/ni
            zi <- qnorm(p1i, lower.tail = FALSE)
            fzi <- dnorm(zi)
            yi <- sqrt(p1i * p2i)/fzi * yi
            yi.t <- ifelse(abs(yi) > 1, sign(yi), yi)
            vi <- 1/(ni - 1) * (p1i * p2i/fzi^2 - (3/2 + (1 - 
                p1i * zi/fzi) * (1 + p2i * zi/fzi)) * yi.t^2 + 
                yi.t^4)
        }
        if (is.element(measure, c("D2OR", "D2ORL"))) {
            yi <- pi/sqrt(3) * di
            vi <- pi^2/3 * (1/n1i + 1/n2i + di^2/(2 * ni))
        }
        if (measure == "D2ORN") {
            yi <- 1.65 * di
            vi <- 1.65^2 * (1/n1i + 1/n2i + di^2/(2 * ni))
        }
    }
    if (is.element(measure, c("COR", "UCOR", "ZCOR"))) {
        mf.ri <- mf[[match("ri", names(mf))]]
        mf.ni <- mf[[match("ni", names(mf))]]
        ri <- eval(mf.ri, data, enclos = sys.frame(sys.parent()))
        ni <- eval(mf.ni, data, enclos = sys.frame(sys.parent()))
        if (!is.null(subset)) {
            ri <- ri[subset]
            ni <- ni[subset]
        }
        if (length(ri) == 0L || length(ni) == 0L) 
            stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments.")
        if (length(ri) != length(ni)) 
            stop("Supplied data vectors are not of the same length.")
        if (any(abs(ri) > 1, na.rm = TRUE)) 
            stop("One or more correlations are > 1 or < -1.")
        if (any(ni < 0, na.rm = TRUE)) 
            stop("One or more sample sizes are negative.")
        if (any(abs(ni) <= 4, na.rm = TRUE)) 
            warning("Cannot estimate sampling variance when ni <= 4.")
        ni.u <- ni
        k <- length(ri)
        if (measure == "COR") {
            yi <- ri
        }
        if (measure == "UCOR") {
            yi <- ri + ri * (1 - ri^2)/(2 * (ni - 4))
            yi[ni <= 4] <- NA
        }
        if (is.element(measure, c("COR", "UCOR"))) {
            if (length(vtype) == 1L) 
                vtype <- rep(vtype, k)
            vi <- rep(NA_real_, k)
            mnwyi <- sum(ni * yi, na.rm = TRUE)/sum(ni, na.rm = TRUE)
            for (i in seq_len(k)) {
                if (vtype[i] == "UB") {
                  vi[i] <- yi[i]^2 - 1 + (ni[i] - 3)/(ni[i] - 
                    2) * ((1 - ri[i]^2) + 2 * (1 - ri[i]^2)^2/ni[i] + 
                    8 * (1 - ri[i]^2)^3/(ni[i] * (ni[i] + 2)) + 
                    48 * (1 - ri[i]^2)^4/(ni[i] * (ni[i] + 2) * 
                      (ni[i] + 4)))
                }
                if (vtype[i] == "LS") {
                  vi[i] <- (1 - yi[i]^2)^2/(ni[i] - 1)
                }
                if (vtype[i] == "HO") 
                  vi[i] <- (1 - mnwyi^2)^2/(ni[i] - 1)
            }
        }
        if (measure == "ZCOR") {
            yi <- 1/2 * log((1 + ri)/(1 - ri))
            vi <- 1/(ni - 3)
        }
        vi[ni <= 4] <- NA
    }
    if (is.element(measure, c("PR", "PLN", "PLO", "PAS", "PFT"))) {
        mf.xi <- mf[[match("xi", names(mf))]]
        mf.mi <- mf[[match("mi", names(mf))]]
        mf.ni <- mf[[match("ni", names(mf))]]
        xi <- eval(mf.xi, data, enclos = sys.frame(sys.parent()))
        mi <- eval(mf.mi, data, enclos = sys.frame(sys.parent()))
        ni <- eval(mf.ni, data, enclos = sys.frame(sys.parent()))
        if (is.null(mi)) 
            mi <- ni - xi
        if (!is.null(subset)) {
            xi <- xi[subset]
            mi <- mi[subset]
        }
        if (length(xi) == 0L || length(mi) == 0L) 
            stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments.")
        if (length(xi) != length(mi)) 
            stop("Supplied data vectors are not all of the same length.")
        if (any(c(xi, mi) < 0, na.rm = TRUE)) 
            stop("One or more counts are negative.")
        ni.u <- xi + mi
        k <- length(xi)
        if (to == "all") {
            xi <- xi + add
            mi <- mi + add
        }
        if (to == "only0") {
            id0 <- c(xi == 0L | mi == 0L)
            id0[is.na(id0)] <- FALSE
            xi[id0] <- xi[id0] + add
            mi[id0] <- mi[id0] + add
        }
        if (to == "if0all") {
            id0 <- c(xi == 0L | mi == 0L)
            id0[is.na(id0)] <- FALSE
            if (any(id0)) {
                xi <- xi + add
                mi <- mi + add
            }
        }
        ni <- xi + mi
        pri <- xi/ni
        if (measure == "PR") {
            yi <- pri
            if (length(vtype) == 1L) 
                vtype <- rep(vtype, k)
            vi <- rep(NA_real_, k)
            mnwpri <- sum(xi, na.rm = TRUE)/sum(ni, na.rm = TRUE)
            for (i in seq_len(k)) {
                if (vtype[i] == "UB") 
                  vi[i] <- pri[i] * (1 - pri[i])/(ni[i] - 1)
                if (vtype[i] == "LS") 
                  vi[i] <- pri[i] * (1 - pri[i])/ni[i]
                if (vtype[i] == "HO") 
                  vi[i] <- mnwpri * (1 - mnwpri)/ni[i]
            }
        }
        if (measure == "PLN") {
            yi <- log(pri)
            if (length(vtype) == 1L) 
                vtype <- rep(vtype, k)
            vi <- rep(NA_real_, k)
            mnwpri <- sum(xi, na.rm = TRUE)/sum(ni, na.rm = TRUE)
            for (i in seq_len(k)) {
                if (vtype[i] == "LS") 
                  vi[i] <- 1/xi[i] - 1/ni[i]
                if (vtype[i] == "HO") 
                  vi[i] <- 1/(mnwpri * ni[i]) - 1/ni[i]
            }
        }
        if (measure == "PLO") {
            yi <- log(pri/(1 - pri))
            if (length(vtype) == 1L) 
                vtype <- rep(vtype, k)
            vi <- rep(NA_real_, k)
            mnwpri <- sum(xi, na.rm = TRUE)/sum(ni, na.rm = TRUE)
            for (i in seq_len(k)) {
                if (vtype[i] == "LS") 
                  vi[i] <- 1/xi[i] + 1/mi[i]
                if (vtype[i] == "HO") 
                  vi[i] <- 1/(mnwpri * ni[i]) + 1/((1 - mnwpri) * 
                    ni[i])
            }
        }
        if (measure == "PAS") {
            yi <- asin(sqrt(pri))
            vi <- 1/(4 * ni)
        }
        if (measure == "PFT") {
            yi <- 1/2 * (asin(sqrt(xi/(ni + 1))) + asin(sqrt((xi + 
                1)/(ni + 1))))
            vi <- 1/(4 * ni + 2)
        }
    }
    if (is.element(measure, c("IR", "IRLN", "IRS", "IRFT"))) {
        mf.xi <- mf[[match("xi", names(mf))]]
        mf.ti <- mf[[match("ti", names(mf))]]
        xi <- eval(mf.xi, data, enclos = sys.frame(sys.parent()))
        ti <- eval(mf.ti, data, enclos = sys.frame(sys.parent()))
        if (!is.null(subset)) {
            xi <- xi[subset]
            ti <- ti[subset]
        }
        if (length(xi) == 0L || length(ti) == 0L) 
            stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments.")
        if (length(xi) != length(ti)) 
            stop("Supplied data vectors are not all of the same length.")
        if (any(xi < 0, na.rm = TRUE)) 
            stop("One or more counts are negative.")
        if (any(ti < 0, na.rm = TRUE)) 
            stop("One or more person-times are negative.")
        ni.u <- ti
        if (to == "all") {
            xi <- xi + add
        }
        if (to == "only0") {
            id0 <- c(xi == 0L)
            id0[is.na(id0)] <- FALSE
            xi[id0] <- xi[id0] + add
        }
        if (to == "if0all") {
            id0 <- c(xi == 0L)
            id0[is.na(id0)] <- FALSE
            if (any(id0)) {
                xi <- xi + add
            }
        }
        iri <- xi/ti
        if (measure == "IR") {
            yi <- iri
            vi <- iri/ti
        }
        if (measure == "IRLN") {
            yi <- log(iri)
            vi <- 1/xi
        }
        if (measure == "IRS") {
            yi <- sqrt(iri)
            vi <- 1/(4 * ti)
        }
        if (measure == "IRFT") {
            yi <- 1/2 * (sqrt(iri) + sqrt(iri + 1/ti))
            vi <- 1/(4 * ti)
        }
    }
    if (is.element(measure, c("MN"))) {
        mf.mi <- mf[[match("mi", names(mf))]]
        mf.sdi <- mf[[match("sdi", names(mf))]]
        mf.ni <- mf[[match("ni", names(mf))]]
        mi <- eval(mf.mi, data, enclos = sys.frame(sys.parent()))
        sdi <- eval(mf.sdi, data, enclos = sys.frame(sys.parent()))
        ni <- eval(mf.ni, data, enclos = sys.frame(sys.parent()))
        if (!is.null(subset)) {
            mi <- mi[subset]
            sdi <- sdi[subset]
            ni <- ni[subset]
        }
        if (length(mi) == 0L || length(sdi) == 0L || length(ni) == 
            0L) 
            stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments.")
        if (!all(length(mi) == c(length(mi), length(sdi), length(ni)))) 
            stop("Supplied data vectors are not all of the same length.")
        if (any(sdi < 0, na.rm = TRUE)) 
            stop("One or more standard deviations are negative.")
        if (any(ni < 0, na.rm = TRUE)) 
            stop("One or more sample sizes are negative.")
        ni.u <- ni
        if (measure == "MN") {
            yi <- mi
            vi <- sdi^2/ni
        }
    }
    if (is.element(measure, c("MC", "SMCC", "SMCR", "SMCRH", 
        "ROMC"))) {
        mf.m1i <- mf[[match("m1i", names(mf))]]
        mf.m2i <- mf[[match("m2i", names(mf))]]
        mf.sd1i <- mf[[match("sd1i", names(mf))]]
        mf.sd2i <- mf[[match("sd2i", names(mf))]]
        mf.ni <- mf[[match("ni", names(mf))]]
        mf.ri <- mf[[match("ri", names(mf))]]
        m1i <- eval(mf.m1i, data, enclos = sys.frame(sys.parent()))
        m2i <- eval(mf.m2i, data, enclos = sys.frame(sys.parent()))
        sd1i <- eval(mf.sd1i, data, enclos = sys.frame(sys.parent()))
        sd2i <- eval(mf.sd2i, data, enclos = sys.frame(sys.parent()))
        ni <- eval(mf.ni, data, enclos = sys.frame(sys.parent()))
        ri <- eval(mf.ri, data, enclos = sys.frame(sys.parent()))
        if (!is.null(subset)) {
            m1i <- m1i[subset]
            m2i <- m2i[subset]
            sd1i <- sd1i[subset]
            sd2i <- sd2i[subset]
            ni <- ni[subset]
            ri <- ri[subset]
        }
        if (is.element(measure, c("MC", "SMCC", "SMCRH", "ROMC"))) {
            if (length(m1i) == 0L || length(m2i) == 0L || length(sd1i) == 
                0L || length(sd2i) == 0L || length(ni) == 0L || 
                length(ri) == 0L) 
                stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments.")
            if (!all(length(m1i) == c(length(m1i), length(m2i), 
                length(sd1i), length(sd2i), length(ni), length(ri)))) 
                stop("Supplied data vectors are not all of the same length.")
            if (any(c(sd1i, sd2i) < 0, na.rm = TRUE)) 
                stop("One or more standard deviations are negative.")
        }
        else {
            if (length(m1i) == 0L || length(m2i) == 0L || length(sd1i) == 
                0L || length(ni) == 0L || length(ri) == 0L) 
                stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments.")
            if (!all(length(m1i) == c(length(m1i), length(m2i), 
                length(sd1i), length(ni), length(ri)))) 
                stop("Supplied data vectors are not all of the same length.")
            if (any(sd1i < 0, na.rm = TRUE)) 
                stop("One or more standard deviations are negative.")
        }
        if (any(abs(ri) > 1, na.rm = TRUE)) 
            stop("One or more correlations are > 1 or < -1.")
        if (any(ni < 0, na.rm = TRUE)) 
            stop("One or more sample sizes are negative.")
        ni.u <- ni
        ni <- ni.u
        mi <- ni - 1
        if (measure == "MC") {
            yi <- m1i - m2i
            vi <- (sd1i^2 + sd2i^2 - 2 * ri * sd1i * sd2i)/ni
        }
        if (measure == "SMCC") {
            warn.before <- getOption("warn")
            options(warn = -1)
            cmi <- .cmicalc(mi)
            options(warn = warn.before)
            sddi <- sqrt(sd1i^2 + sd2i^2 - 2 * ri * sd1i * sd2i)
            yi <- cmi * (m1i - m2i)/sddi
            vi <- 1/ni + yi^2/(2 * ni)
        }
        if (measure == "SMCR") {
            warn.before <- getOption("warn")
            options(warn = -1)
            cmi <- .cmicalc(mi)
            options(warn = warn.before)
            yi <- cmi * (m1i - m2i)/sd1i
            vi <- 2 * (1 - ri)/ni + yi^2/(2 * ni)
        }
        if (measure == "SMCRH") {
            warn.before <- getOption("warn")
            options(warn = -1)
            cmi <- .cmicalc(mi)
            options(warn = warn.before)
            vardi <- sd1i^2 + sd2i^2 - 2 * ri * sd1i * sd2i
            yi <- cmi * (m1i - m2i)/sd1i
            vi <- vardi/(sd1i^2 * (ni - 1)) + yi^2/(2 * (ni - 
                1))
            vi <- cmi^2 * vi
        }
        if (measure == "ROMC") {
            yi <- log(m1i/m2i)
            vi <- sd1i^2/(ni * m1i^2) + sd2i^2/(ni * m2i^2) - 
                2 * ri * sd1i * sd2i/(m1i * m2i * ni)
        }
    }
    if (is.element(measure, c("ARAW", "AHW", "ABT"))) {
        mf.ai <- mf[[match("ai", names(mf))]]
        mf.mi <- mf[[match("mi", names(mf))]]
        mf.ni <- mf[[match("ni", names(mf))]]
        ai <- eval(mf.ai, data, enclos = sys.frame(sys.parent()))
        mi <- eval(mf.mi, data, enclos = sys.frame(sys.parent()))
        ni <- eval(mf.ni, data, enclos = sys.frame(sys.parent()))
        if (!is.null(subset)) {
            ai <- ai[subset]
            mi <- mi[subset]
            ni <- ni[subset]
        }
        if (length(ai) == 0L || length(mi) == 0L || length(ni) == 
            0L) 
            stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments.")
        if (!all(length(ai) == c(length(ai), length(mi), length(ni)))) 
            stop("Supplied data vectors are not all of the same length.")
        if (any(ai > 1, na.rm = TRUE)) 
            stop("One or more alpha values are > 1.")
        if (any(mi < 2, na.rm = TRUE)) 
            stop("One or more mi values are < 2.")
        if (any(ni < 0, na.rm = TRUE)) 
            stop("One or more sample sizes are negative.")
        ni.u <- ni
        if (measure == "ARAW") {
            yi <- ai
            vi <- 2 * mi * (1 - ai)^2/((mi - 1) * (ni - 2))
        }
        if (measure == "AHW") {
            yi <- 1 - (1 - ai)^(1/3)
            vi <- 18 * mi * (ni - 1) * (1 - ai)^(2/3)/((mi - 
                1) * (9 * ni - 11)^2)
        }
        if (measure == "ABT") {
            yi <- -log(1 - ai)
            vi <- 2 * mi/((mi - 1) * (ni - 2))
        }
    }
    is.inf <- is.infinite(yi) | is.infinite(vi)
    if (any(is.inf)) {
        warning("Some yi and/or vi values equal to +-Inf. Recoded to NAs.")
        yi[is.inf] <- NA
        vi[is.inf] <- NA
    }
    is.NaN <- is.nan(yi) | is.nan(vi)
    if (any(is.NaN)) {
        yi[is.NaN] <- NA
        vi[is.NaN] <- NA
    }
    vi[vi < 0] <- NA
    if (!is.null(slab)) {
        if (!is.null(subset)) 
            slab <- slab[subset]
        if (anyNA(slab)) 
            stop("NAs in study labels.")
        if (anyDuplicated(slab)) 
            slab <- make.unique(as.character(slab))
        if (length(slab) != length(yi)) 
            stop("Study labels not of same length as data.")
        attr(yi, "slab") <- slab
    }
    if (!is.null(subset)) {
        if (!no.data) 
            data <- data[subset, , drop = FALSE]
    }
    attr(yi, "measure") <- measure
    if (!no.data && append) {
        dat <- data.frame(data)
        if (replace) {
            dat[[var.names[1]]] <- yi
            dat[[var.names[2]]] <- vi
            if (add.measure) {
                dat[[var.names[3]]] <- ""
                dat[[var.names[3]]][!is.na(yi)] <- measure
            }
            attr(dat[[var.names[1]]], "ni") <- ni.u
        }
        else {
            if (is.element(var.names[1], names(dat))) {
                is.na.yi <- is.na(dat[[var.names[1]]])
                dat[[var.names[1]]][is.na.yi] <- yi[is.na.yi]
                attributes(dat[[var.names[1]]])$ni[is.na.yi] <- ni.u[is.na.yi]
            }
            else {
                dat[[var.names[1]]] <- yi
                attr(dat[[var.names[1]]], "ni") <- ni.u
            }
            if (is.element(var.names[2], names(dat))) {
                is.na.vi <- is.na(dat[[var.names[2]]])
                dat[[var.names[2]]][is.na.vi] <- vi[is.na.vi]
            }
            else {
                dat[[var.names[2]]] <- vi
            }
            if (add.measure) {
                if (is.element(var.names[3], names(dat))) {
                  is.na.measure <- c(dat[[var.names[3]]] == "") & 
                    !is.na(yi)
                  dat[[var.names[3]]][is.na.measure] <- measure
                }
                else {
                  dat[[var.names[3]]] <- ""
                  dat[[var.names[3]]][!is.na(yi)] <- measure
                }
            }
        }
    }
    else {
        if (add.measure) {
            dat <- data.frame(yi, vi)
            dat$measure <- ""
            dat$measure[!is.na(yi)] <- measure
            names(dat) <- var.names
        }
        else {
            dat <- data.frame(yi, vi)
            names(dat) <- var.names[1:2]
        }
        attr(dat[, 1], "ni") <- ni.u
    }
    attr(dat, "digits") <- digits
    attr(dat, "yi.names") <- unique(c(var.names[1], attr(data, 
        "yi.names")))
    attr(dat, "vi.names") <- unique(c(var.names[2], attr(data, 
        "vi.names")))
    attr(dat, "sei.names") <- attr(data, "sei.names")
    attr(dat, "zi.names") <- attr(data, "zi.names")
    attr(dat, "ci.lb.names") <- attr(data, "ci.lb.names")
    attr(dat, "ci.ub.names") <- attr(data, "ci.ub.names")
    class(dat) <- c("escalc", "data.frame")
    return(dat)
}
