rma.mv <-
function (yi, V, W, mods, random, struct = "CS", intercept = TRUE, 
    data, slab, subset, method = "REML", tdist = FALSE, level = 95, 
    digits = 4, btt, R, Rscale = "cor", sigma2, tau2, rho, gamma2, 
    phi, sparse = FALSE, verbose = FALSE, control) 
{
    if (!is.element(method, c("FE", "ML", "REML"))) 
        stop("Unknown 'method' specified.")
    if (any(!is.element(struct, c("CS", "HCS", "UN", "AR", "HAR", 
        "UNHO", "ID", "DIAG")))) 
        stop("Unknown 'struct' specified.")
    if (length(struct) == 1) 
        struct <- c(struct, struct)
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    if (missing(random)) 
        random <- NULL
    if (missing(btt)) 
        btt <- NULL
    if (missing(R)) 
        R <- NULL
    if (missing(sigma2)) 
        sigma2 <- NULL
    if (missing(tau2)) 
        tau2 <- NULL
    if (missing(rho)) 
        rho <- NULL
    if (missing(gamma2)) 
        gamma2 <- NULL
    if (missing(phi)) 
        phi <- NULL
    if (missing(control)) 
        control <- list()
    knha <- tdist
    if (is.character(Rscale)) 
        Rscale <- match.arg(Rscale, c("none", "cor", "cor0", 
            "cov0"))
    if (is.logical(Rscale)) 
        Rscale <- ifelse(Rscale, "cor", "none")
    if (is.numeric(Rscale)) {
        Rscale <- round(Rscale)
        if (Rscale > 3 | Rscale < 0) 
            stop("Unknown 'Rscale' value specified.")
        Rscale <- switch(as.character(Rscale), `0` = "none", 
            `1` = "cor", `2` = "cor0", `3` = "cov0")
    }
    very.verbose <- ifelse(!is.logical(verbose) && verbose > 
        1, TRUE, FALSE)
    if (very.verbose) 
        message("Extracting yi/V values ...")
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
    mf.V <- mf[[match("V", names(mf))]]
    mf.W <- mf[[match("W", names(mf))]]
    mf.ni <- mf[[match("ni", names(mf))]]
    mf.slab <- mf[[match("slab", names(mf))]]
    mf.subset <- mf[[match("subset", names(mf))]]
    mf.mods <- mf[[match("mods", names(mf))]]
    yi <- eval(mf.yi, data, enclos = sys.frame(sys.parent()))
    V <- eval(mf.V, data, enclos = sys.frame(sys.parent()))
    W <- eval(mf.W, data, enclos = sys.frame(sys.parent()))
    ni <- eval(mf.ni, data, enclos = sys.frame(sys.parent()))
    slab <- eval(mf.slab, data, enclos = sys.frame(sys.parent()))
    subset <- eval(mf.subset, data, enclos = sys.frame(sys.parent()))
    mods <- eval(mf.mods, data, enclos = sys.frame(sys.parent()))
    is.formula <- FALSE
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
    measure <- "GEN"
    if (!is.null(attr(yi, "measure"))) 
        measure <- attr(yi, "measure")
    attr(yi, "measure") <- measure
    if (is.null(V)) 
        stop("Need to specify V argument.")
    if (is.list(V)) {
        rows <- sapply(V, NROW)
        cols <- sapply(V, NCOL)
        if (any(rows != cols)) 
            stop("List elements in V must be square matrices.")
        if (sparse) {
            V <- bdiag(V)
        }
        else {
            V <- bldiag(V)
        }
    }
    if (is.vector(V) || nrow(V) == 1L || ncol(V) == 1L) 
        V <- diag(as.vector(V), nrow = k, ncol = k)
    if (is.data.frame(V)) 
        V <- as.matrix(V)
    V <- unname(V)
    if (dim(V)[1] != dim(V)[2]) 
        stop("V must be a square matrix.")
    if (!isSymmetric(V)) 
        stop("V must be a symmetric matrix.")
    if (dim(V)[1] != k) 
        stop("Length of yi and length/dimensions of V are not the same.")
    if (sparse && class(V) == "matrix") 
        V <- Matrix(V, sparse = TRUE)
    if (!is.null(W)) {
        if (is.vector(W) || nrow(W) == 1L || ncol(W) == 1L) {
            W <- as.vector(W)
            if (length(W) == 1L) 
                W <- rep(W, k)
            A <- diag(W, nrow = length(W), ncol = length(W))
        }
        else {
            A <- W
        }
        if (is.data.frame(A)) 
            A <- as.matrix(A)
        A <- unname(A)
        if (dim(A)[1] != dim(A)[2]) 
            stop("W must be a square matrix.")
        if (!isSymmetric(A)) 
            stop("W must be a symmetric matrix.")
        if (dim(A)[1] != k) 
            stop("Length of yi and length/dimensions of W are not the same.")
        if (sparse && class(A) == "matrix") 
            A <- Matrix(A, sparse = TRUE)
    }
    else {
        A <- NULL
    }
    if (is.null(ni) && !is.null(attr(yi, "ni"))) 
        ni <- attr(yi, "ni")
    if (!is.null(ni) && length(ni) != k) 
        ni <- NULL
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
    if (method != "FE" && !is.null(random)) {
        if (very.verbose) 
            message("Processing 'random' argument ...")
        if (!is.list(random)) 
            random <- list(random)
        has.slash <- sapply(random, function(f) grepl("/", paste0(f, 
            collapse = "")))
        random.plus <- lapply(random, function(f) formula(sub("\\|", 
            "+", paste0(f, collapse = ""))))
        mf.r <- lapply(random.plus, model.frame, data = data, 
            na.action = na.pass)
        mf.r.ncols <- sapply(mf.r, ncol)
        for (j in 1:length(has.slash)) {
            if (!has.slash[j]) 
                next
            for (p in mf.r.ncols[j]:1) {
                mf.r[[j]][, p] <- interaction(mf.r[[j]][1:p], 
                  drop = TRUE)
                colnames(mf.r[[j]])[p] <- paste(colnames(mf.r[[j]])[1:p], 
                  collapse = "/")
            }
        }
        if (any(has.slash)) {
            if (length(mf.r) == 1) {
                mf.r <- lapply(seq(ncol(mf.r[[1]])), function(x) mf.r[[1]][x])
            }
            else {
                mf.r <- unlist(mapply(function(mf, sl) if (sl) 
                  lapply(seq(mf), function(x) mf[x])
                else list(mf), mf.r, has.slash), recursive = FALSE, 
                  use.names = FALSE)
            }
            mf.r.ncols <- sapply(mf.r, ncol)
        }
        if (any(mf.r.ncols > 2)) 
            stop("No more than two elements allowed in each formula of the 'random' argument.")
        if (sum(mf.r.ncols == 2) > 2) 
            stop("Only up to two formulas with two elements allowed in the 'random' argument.")
        mf.s <- mf.r[which(mf.r.ncols == 1)]
        mf.g <- mf.r[[which(mf.r.ncols == 2)[1]]]
        mf.h <- mf.r[[which(mf.r.ncols == 2)[2]]]
        if (length(mf.s) == 0) 
            mf.s <- NULL
        withS <- !is.null(mf.s)
        withG <- !is.null(mf.g)
        withH <- !is.null(mf.h)
        mf.r.nrows <- sapply(mf.r, nrow)
        if (any(mf.r.nrows != k)) 
            stop("Length of variables specified via the 'random' argument does not match length of the data.")
        mf.r <- lapply(random.plus, get_all_vars, data = data)
    }
    else {
        mf.s <- NULL
        mf.g <- NULL
        mf.h <- NULL
        mf.r <- NULL
        withS <- FALSE
        withG <- FALSE
        withH <- FALSE
    }
    ids <- seq_len(k)
    if (very.verbose) 
        message("Generating/extracting study labels ...")
    if (is.null(slab)) {
        if (!is.null(attr(yi, "slab"))) 
            slab <- attr(yi, "slab")
        if (is.null(slab) && length(slab) != k) 
            slab <- NULL
    }
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
        yi <- yi[subset]
        V <- V[subset, subset, drop = FALSE]
        A <- A[subset, subset, drop = FALSE]
        ni <- ni[subset]
        mods <- mods[subset, , drop = FALSE]
        slab <- slab[subset]
        mf.s <- lapply(mf.s, function(x) x[subset, , drop = FALSE])
        mf.g <- mf.g[subset, , drop = FALSE]
        mf.h <- mf.h[subset, , drop = FALSE]
        mf.r <- lapply(mf.r, function(x) x[subset, , drop = FALSE])
        ids <- ids[subset]
        k <- length(yi)
        attr(yi, "measure") <- measure
        attr(yi, "ni") <- ni
    }
    if (anyDuplicated(slab)) 
        slab <- make.unique(as.character(slab))
    attr(yi, "slab") <- slab
    vi <- diag(V)
    if (any(vi <= 0, na.rm = TRUE)) {
        allvipos <- FALSE
        warning("There are outcomes with non-positive sampling variances.")
        vi.neg <- vi < 0
        if (any(vi.neg, na.rm = TRUE)) {
            V[vi.neg, , drop = FALSE] <- 0
            V[, vi.neg, drop = FALSE] <- 0
            vi[vi.neg] <- 0
            warning("Negative sampling variances constrained to zero.")
        }
    }
    else {
        allvipos <- TRUE
    }
    yi.f <- yi
    vi.f <- vi
    V.f <- V
    W.f <- A
    ni.f <- ni
    mods.f <- mods
    mf.g.f <- mf.g
    mf.h.f <- mf.h
    k.f <- k
    if (withS) {
        if (very.verbose) 
            message(paste0("Processing '", paste(as.character(random[mf.r.ncols == 
                1]), collapse = ", "), "' term(s) ..."))
        s.names <- sapply(mf.s, names)
        mf.s <- lapply(mf.s, function(x) factor(x[[1]]))
        if (any(sapply(lapply(mf.s, is.na), any))) 
            stop("No NAs allowed in variables specified in the 'random' argument.")
        sigma2s <- length(mf.s)
        if (is.null(sigma2)) 
            sigma2 <- rep(NA_real_, sigma2s)
        if (length(sigma2) == 1) 
            sigma2 <- rep(sigma2, sigma2s)
        if (length(sigma2) != sigma2s) 
            stop(paste("Length of 'sigma2' argument (", length(sigma2), 
                ") does not match actual number of variance components (", 
                sigma2s, ").", sep = ""))
        if (any(sigma2 < 0, na.rm = TRUE)) 
            stop("Specified value(s) of 'sigma2' must be non-negative.")
        s.nlevels <- sapply(mf.s, nlevels)
        s.levels <- lapply(mf.s, levels)
        if (is.null(R)) {
            withR <- FALSE
            Rfix <- rep(FALSE, sigma2s)
        }
        else {
            if (very.verbose) 
                message("Processing 'R' argument ...")
            withR <- TRUE
            if (is.data.frame(R) || !is.list(R)) 
                R <- list(R)
            if (is.null(names(R)) || any(nchar(names(R)) == 0)) 
                stop("Argument 'R' must be a *named* list.")
            R <- R[!sapply(R, is.null)]
            R <- lapply(R, as.matrix)
            R <- R[s.names]
            names(R) <- s.names
            Rfix <- !sapply(R, is.null)
            if (any(Rfix)) {
                if (any(sapply(R[Rfix], function(x) dim(x)[1] != 
                  dim(x)[2]))) 
                  stop("Elements of 'R' must be square matrices.")
                if (any(sapply(R[Rfix], function(x) !isSymmetric(unname(x))))) 
                  stop("Elements of 'R' must be symmetric matrices.")
                for (j in 1:length(R)) {
                  if (!Rfix[j]) 
                    next
                  if (is.null(rownames(R[[j]]))) 
                    rownames(R[[j]]) <- colnames(R[[j]])
                  if (is.null(colnames(R[[j]]))) 
                    colnames(R[[j]]) <- rownames(R[[j]])
                  if (is.null(colnames(R[[j]]))) 
                    stop("Elements of 'R' must have dimension names.")
                }
                R[Rfix] <- lapply(R[Rfix], function(x) x[!duplicated(colnames(x)), 
                  !duplicated(colnames(x)), drop = FALSE])
                if (any(sapply(R[Rfix], function(x) length(colnames(x)) != 
                  length(unique(colnames(x)))))) 
                  stop("Each element of 'R' must have unique dimension names.")
                for (j in 1:length(R)) {
                  if (!Rfix[j]) 
                    next
                  if (anyNA(R[[j]])) 
                    stop("No missing values allowed in matrix specified via 'R'.")
                  if (any(!is.element(s.levels[[j]], colnames(R[[j]])))) 
                    stop(paste0("There are levels in '", s.names[j], 
                      "' for which there are no rows/columns in the corresponding 'R' matrix."))
                  if (any(!is.element(colnames(R[[j]]), s.levels[[j]]))) 
                    warning(paste0("There are rows/columns in the 'R' matrix for '", 
                      s.names[j], "' for which there are no data."))
                }
            }
            else {
                warning("Argument 'R' specified, but list name(s) not in 'random'.")
                withR <- FALSE
                Rfix <- rep(FALSE, sigma2s)
                R <- NULL
            }
        }
    }
    else {
        sigma2s <- 1
        sigma2 <- 0
        s.nlevels <- NULL
        s.levels <- NULL
        s.names <- NULL
        withR <- FALSE
        Rfix <- FALSE
        R <- NULL
    }
    if (withG) {
        if (very.verbose) 
            message(paste0("Processing '", as.character(random[mf.r.ncols == 
                2][1]), "' term ..."))
        g.names <- names(mf.g)
        if (is.element(struct[1], c("CS", "HCS", "UN", "ID", 
            "DIAG", "UNHO")) && !is.factor(mf.g.f[[1]]) && !is.character(mf.g.f[[1]])) 
            stop("Inner variable in (~ inner | outer) must be a factor or character variable.")
        mf.g.f <- data.frame(inner = factor(mf.g.f[[1]]), outer = factor(mf.g.f[[2]]))
        mf.g <- data.frame(inner = factor(mf.g[[1]]), outer = factor(mf.g[[2]]))
        if (anyNA(mf.g)) 
            stop("No NAs allowed in variables specified in the 'random' argument.")
        g.nlevels <- c(nlevels(mf.g[[1]]), nlevels(mf.g[[2]]))
        g.levels <- list(levels(mf.g[[1]]), levels(mf.g[[2]]))
        if (struct[1] == "CS" || struct[1] == "ID") {
            tau2s <- 1
            rhos <- 1
        }
        if (struct[1] == "HCS" || struct[1] == "DIAG") {
            tau2s <- g.nlevels[1]
            rhos <- 1
        }
        if (struct[1] == "UN") {
            tau2s <- g.nlevels[1]
            rhos <- ifelse(g.nlevels[1] > 1, g.nlevels[1] * (g.nlevels[1] - 
                1)/2, 1)
        }
        if (struct[1] == "AR") {
            tau2s <- 1
            rhos <- 1
        }
        if (struct[1] == "HAR") {
            tau2s <- g.nlevels[1]
            rhos <- 1
        }
        if (struct[1] == "UNHO") {
            tau2s <- 1
            rhos <- ifelse(g.nlevels[1] > 1, g.nlevels[1] * (g.nlevels[1] - 
                1)/2, 1)
        }
        if (is.null(tau2)) 
            tau2 <- rep(NA_real_, tau2s)
        if (is.null(rho)) 
            rho <- rep(NA_real_, rhos)
        if (length(tau2) == 1 && is.element(struct[1], c("HCS", 
            "UN", "DIAG", "HAR"))) 
            tau2 <- rep(tau2, tau2s)
        if (length(rho) == 1 && is.element(struct[1], c("UN", 
            "UNHO"))) 
            rho <- rep(rho, rhos)
        if (length(tau2) != tau2s) 
            stop(paste("Length of 'tau2' argument (", length(tau2), 
                ") does not match actual number of variance components (", 
                tau2s, ").", sep = ""))
        if (length(rho) != rhos) 
            stop(paste("Length of 'rho' argument (", length(rho), 
                ") does not match actual number of correlations (", 
                rhos, ").", sep = ""))
        if (any(tau2 < 0, na.rm = TRUE)) 
            stop("Specified value(s) of 'tau2' must be non-negative.")
        if (any(rho > 1 | rho < -1, na.rm = TRUE)) 
            stop("Specified value(s) of 'rho' must be in [-1,1].")
        if (g.nlevels[1] == 1) {
            Z.G1 <- cbind(rep(1, k))
        }
        else {
            if (sparse) {
                Z.G1 <- sparse.model.matrix(~mf.g[[1]] - 1)
            }
            else {
                Z.G1 <- model.matrix(~mf.g[[1]] - 1)
            }
        }
        if (g.nlevels[2] == 1) {
            Z.G2 <- cbind(rep(1, k))
        }
        else {
            if (sparse) {
                Z.G2 <- sparse.model.matrix(~mf.g[[2]] - 1)
            }
            else {
                Z.G2 <- model.matrix(~mf.g[[2]] - 1)
            }
        }
    }
    else {
        tau2s <- 1
        rhos <- 1
        tau2 <- 0
        rho <- 0
        Z.G1 <- NULL
        Z.G2 <- NULL
        g.nlevels <- NULL
        g.levels <- NULL
        g.names <- NULL
    }
    if (withH) {
        if (very.verbose) 
            message(paste0("Processing '", as.character(random[mf.r.ncols == 
                2][2]), "' term ..."))
        h.names <- names(mf.h)
        if (is.element(struct[2], c("CS", "HCS", "UN", "ID", 
            "DIAG", "UNHO")) && !is.factor(mf.h.f[[1]]) && !is.character(mf.h.f[[1]])) 
            stop("Inner variable in (~ inner | outer) must be a factor or character variable.")
        mf.h.f <- data.frame(inner = factor(mf.h.f[[1]]), outer = factor(mf.h.f[[2]]))
        mf.h <- data.frame(inner = factor(mf.h[[1]]), outer = factor(mf.h[[2]]))
        if (anyNA(mf.h)) 
            stop("No NAs allowed in variables specified in the 'random' argument.")
        h.nlevels <- c(nlevels(mf.h[[1]]), nlevels(mf.h[[2]]))
        h.levels <- list(levels(mf.h[[1]]), levels(mf.h[[2]]))
        if (struct[2] == "CS" || struct[2] == "ID") {
            gamma2s <- 1
            phis <- 1
        }
        if (struct[2] == "HCS" || struct[2] == "DIAG") {
            gamma2s <- h.nlevels[1]
            phis <- 1
        }
        if (struct[2] == "UN") {
            gamma2s <- h.nlevels[1]
            phis <- ifelse(h.nlevels[1] > 1, h.nlevels[1] * (h.nlevels[1] - 
                1)/2, 1)
        }
        if (struct[2] == "AR") {
            gamma2s <- 1
            phis <- 1
        }
        if (struct[2] == "HAR") {
            gamma2s <- h.nlevels[1]
            phis <- 1
        }
        if (struct[2] == "UNHO") {
            gamma2s <- 1
            phis <- ifelse(h.nlevels[1] > 1, h.nlevels[1] * (h.nlevels[1] - 
                1)/2, 1)
        }
        if (is.null(gamma2)) 
            gamma2 <- rep(NA_real_, gamma2s)
        if (is.null(phi)) 
            phi <- rep(NA_real_, phis)
        if (length(gamma2) == 1 && is.element(struct[2], c("HCS", 
            "UN", "DIAG", "HAR"))) 
            gamma2 <- rep(gamma2, gamma2s)
        if (length(phi) == 1 && is.element(struct[2], c("UN", 
            "UNHO"))) 
            phi <- rep(phi, phis)
        if (length(gamma2) != gamma2s) 
            stop(paste("Length of 'gamma2' argument (", length(gamma2), 
                ") does not match actual number of variance components (", 
                gamma2s, ").", sep = ""))
        if (length(phi) != phis) 
            stop(paste("Length of 'phi' argument (", length(phi), 
                ") does not match actual number of correlations (", 
                phis, ").", sep = ""))
        if (any(gamma2 < 0, na.rm = TRUE)) 
            stop("Specified value(s) of 'gamma2' must be non-negative.")
        if (any(phi > 1 | phi < -1, na.rm = TRUE)) 
            stop("Specified value(s) of 'phi' must be in [-1,1].")
        if (h.nlevels[1] == 1) {
            Z.H1 <- cbind(rep(1, k))
        }
        else {
            if (sparse) {
                Z.H1 <- sparse.model.matrix(~mf.h[[1]] - 1)
            }
            else {
                Z.H1 <- model.matrix(~mf.h[[1]] - 1)
            }
        }
        if (h.nlevels[2] == 1) {
            Z.H2 <- cbind(rep(1, k))
        }
        else {
            if (sparse) {
                Z.H2 <- sparse.model.matrix(~mf.h[[2]] - 1)
            }
            else {
                Z.H2 <- model.matrix(~mf.h[[2]] - 1)
            }
        }
    }
    else {
        gamma2s <- 1
        phis <- 1
        gamma2 <- 0
        phi <- 0
        Z.H1 <- NULL
        Z.H2 <- NULL
        h.nlevels <- NULL
        h.levels <- NULL
        h.names <- NULL
    }
    Vlt <- V
    Vlt[upper.tri(Vlt)] <- 0
    has.na <- is.na(yi) | if (is.null(mods)) 
        FALSE
    else apply(is.na(mods), 1, any) | apply(is.na(Vlt), 1, any) | 
        if (is.null(A)) 
            FALSE
        else apply(is.na(A), 1, any)
    if (any(has.na)) {
        if (very.verbose) 
            message("Handling NAs ...")
        not.na <- !has.na
        if (na.act == "na.omit" || na.act == "na.exclude" || 
            na.act == "na.pass") {
            yi <- yi[not.na]
            V <- V[not.na, not.na, drop = FALSE]
            A <- A[not.na, not.na, drop = FALSE]
            vi <- vi[not.na]
            ni <- ni[not.na]
            mods <- mods[not.na, , drop = FALSE]
            mf.s <- lapply(mf.s, function(x) x[not.na])
            mf.g <- mf.g[not.na, , drop = FALSE]
            mf.h <- mf.h[not.na, , drop = FALSE]
            mf.r <- lapply(mf.r, function(x) x[not.na, , drop = FALSE])
            Z.G1 <- Z.G1[not.na, , drop = FALSE]
            Z.G2 <- Z.G2[not.na, , drop = FALSE]
            Z.H1 <- Z.H1[not.na, , drop = FALSE]
            Z.H2 <- Z.H2[not.na, , drop = FALSE]
            k <- length(yi)
            warning("Rows with NAs omitted from model fitting.")
            attr(yi, "measure") <- measure
            attr(yi, "ni") <- ni
        }
        if (na.act == "na.fail") 
            stop("Missing values in data.")
    }
    else {
        not.na <- rep(TRUE, k)
    }
    if (k <= 1) 
        stop("Processing terminated since k <= 1.")
    if (any(eigen(V, symmetric = TRUE, only.values = TRUE)$values <= 
        .Machine$double.eps)) 
        warning("V appears to be not positive definite.")
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
    if (withS) {
        mf.s <- lapply(mf.s, factor)
        s.nlevels <- sapply(mf.s, nlevels)
        s.levels <- lapply(mf.s, levels)
        if (any(is.na(sigma2) & s.nlevels == 1)) {
            sigma2[is.na(sigma2) & s.nlevels == 1] <- 0
            warning("Single-level factor(s) found in 'random' argument. Corresponding sigma2 value(s) fixed to 0.")
        }
        Z.S <- vector(mode = "list", length = sigma2s)
        for (j in 1:sigma2s) {
            if (s.nlevels[j] == 1) {
                Z.S[[j]] <- cbind(rep(1, k))
            }
            else {
                if (sparse) {
                  Z.S[[j]] <- sparse.model.matrix(~mf.s[[j]] - 
                    1)
                }
                else {
                  Z.S[[j]] <- model.matrix(~mf.s[[j]] - 1)
                }
            }
        }
    }
    else {
        Z.S <- NULL
    }
    if (withR) {
        for (j in 1:length(R)) {
            if (!Rfix[j]) 
                next
            R[[j]] <- R[[j]][s.levels[[j]], s.levels[[j]]]
        }
        if (Rscale == "cor" || Rscale == "cor0") 
            R[Rfix] <- lapply(R[Rfix], cov2cor)
        if (Rscale == "cor0") 
            R[Rfix] <- lapply(R[Rfix], function(x) (x - min(x))/(1 - 
                min(x)))
        if (Rscale == "cov0") 
            R[Rfix] <- lapply(R[Rfix], function(x) (x - min(x)))
    }
    if (withS) {
        D.S <- vector(mode = "list", length = sigma2s)
        for (j in seq_len(sigma2s)) {
            if (Rfix[j]) {
                if (sparse) {
                  D.S[[j]] <- Z.S[[j]] %*% Matrix(R[[j]], sparse = TRUE) %*% 
                    t(Z.S[[j]])
                }
                else {
                  D.S[[j]] <- Z.S[[j]] %*% R[[j]] %*% t(Z.S[[j]])
                }
            }
            else {
                D.S[[j]] <- tcrossprod(Z.S[[j]])
            }
        }
    }
    else {
        D.S <- NULL
    }
    if (withG) {
        g.nlevels.f <- g.nlevels
        g.levels.f <- g.levels
        mf.g <- data.frame(inner = factor(mf.g[[1]]), outer = factor(mf.g[[2]]))
        g.nlevels <- c(nlevels(mf.g[[1]]), nlevels(mf.g[[2]]))
        g.levels <- list(levels(mf.g[[1]]), levels(mf.g[[2]]))
        g.levels.r <- !is.element(g.levels.f[[1]], g.levels[[1]])
        if (any(g.levels.r)) 
            warning("One or more levels of inner factor removed due to NAs.")
        if (is.element(struct[1], c("ID", "DIAG"))) 
            rho <- 0
        if (g.nlevels[1] == 1 && is.element(struct[1], c("CS", 
            "HCS", "AR", "HAR")) && is.na(rho)) {
            rho <- 0
            warning("Inner factor has only a single level, so fixed value of 'rho' to 0.")
        }
        g.levels.k <- table(factor(mf.g[[1]], levels = g.levels.f[[1]]))
        if (is.element(struct[1], c("HCS", "UN", "DIAG", "HAR"))) {
            if (any(is.na(tau2) & g.levels.k == 1)) {
                tau2[is.na(tau2) & g.levels.k == 1] <- 0
                warning("Inner factor has k=1 for one or more levels. Corresponding tau2 value(s) fixed to 0.")
            }
        }
        g.levels.comb.k <- crossprod(Z.G2, Z.G1)
        g.levels.comb.k <- split(g.levels.comb.k, 1:nrow(g.levels.comb.k))
        if (all(unlist(lapply(g.levels.comb.k, sum)) == 1)) {
            if (is.element(struct[1], c("CS", "HCS", "AR", "HAR")) && 
                is.na(rho)) {
                rho <- 0
                warning("Each level of the outer factor contains only a single level of the inner factor, so fixed value of 'rho' to 0.")
            }
        }
        g.levels.comb.k <- lapply(g.levels.comb.k, function(x) outer(x, 
            x, FUN = "&"))
        g.levels.comb.k <- lapply(g.levels.comb.k, function(x) ifelse(x, 
            1, 0))
        g.levels.comb.k <- Reduce("+", g.levels.comb.k)
        g.levels.comb.k <- g.levels.comb.k[upper.tri(g.levels.comb.k)]
        if (is.element(struct[1], c("UN", "UNHO")) && any(g.levels.comb.k == 
            0 & is.na(rho))) {
            rho[g.levels.comb.k == 0] <- 0
            warning("Some combinations of the levels of the inner factor never occurred. Corresponding 'rho' value(s) fixed to 0.")
        }
        if (is.element(struct[1], c("UN", "UNHO")) && g.nlevels.f[1] == 
            1 && is.na(rho)) {
            rho <- 0
            warning("Inner factor has only a single level, so fixed value of 'rho' to 0.")
        }
        if (struct[1] == "CS") {
            G <- matrix(rho * tau2, nrow = g.nlevels.f[1], ncol = g.nlevels.f[1])
            diag(G) <- tau2
        }
        if (struct[1] == "HCS") {
            G <- matrix(rho, nrow = g.nlevels.f[1], ncol = g.nlevels.f[1])
            diag(G) <- 1
            G <- diag(sqrt(tau2), nrow = g.nlevels.f[1], ncol = g.nlevels.f[1]) %*% 
                G %*% diag(sqrt(tau2), nrow = g.nlevels.f[1], 
                ncol = g.nlevels.f[1])
            diag(G) <- tau2
        }
        if (struct[1] == "UN") {
            G <- .con.vcov.UN(tau2, rho)
        }
        if (struct[1] == "ID" || struct[1] == "DIAG") {
            G <- diag(tau2, nrow = g.nlevels.f[1], ncol = g.nlevels.f[1])
        }
        if (struct[1] == "UNHO") {
            G <- matrix(NA_real_, nrow = g.nlevels.f[1], ncol = g.nlevels.f[1])
            G[upper.tri(G)] <- rho
            G[lower.tri(G)] <- t(G)[lower.tri(G)]
            diag(G) <- 1
            G <- diag(sqrt(rep(tau2, g.nlevels.f[1])), nrow = g.nlevels.f[1], 
                ncol = g.nlevels.f[1]) %*% G %*% diag(sqrt(rep(tau2, 
                g.nlevels.f[1])), nrow = g.nlevels.f[1], ncol = g.nlevels.f[1])
            diag(G) <- tau2
        }
        if (struct[1] == "AR") {
            if (is.na(rho)) {
                G <- matrix(NA_real_, nrow = g.nlevels.f[1], 
                  ncol = g.nlevels.f[1])
            }
            else {
                if (g.nlevels.f[1] > 1) {
                  G <- toeplitz(ARMAacf(ar = rho, lag.max = g.nlevels.f[1] - 
                    1))
                }
                else {
                  G <- diag(1)
                }
            }
            G <- diag(sqrt(rep(tau2, g.nlevels.f[1])), nrow = g.nlevels.f[1], 
                ncol = g.nlevels.f[1]) %*% G %*% diag(sqrt(rep(tau2, 
                g.nlevels.f[1])), nrow = g.nlevels.f[1], ncol = g.nlevels.f[1])
            diag(G) <- tau2
        }
        if (struct[1] == "HAR") {
            if (is.na(rho)) {
                G <- matrix(NA_real_, nrow = g.nlevels.f[1], 
                  ncol = g.nlevels.f[1])
            }
            else {
                if (g.nlevels.f[1] > 1) {
                  G <- toeplitz(ARMAacf(ar = rho, lag.max = g.nlevels.f[1] - 
                    1))
                }
                else {
                  G <- diag(1)
                }
            }
            G <- diag(sqrt(tau2), nrow = g.nlevels.f[1], ncol = g.nlevels.f[1]) %*% 
                G %*% diag(sqrt(tau2), nrow = g.nlevels.f[1], 
                ncol = g.nlevels.f[1])
            diag(G) <- tau2
        }
        if (any(g.levels.r) && is.element(struct[1], c("CS", 
            "AR", "ID"))) {
            G[g.levels.r, ] <- 0
            G[, g.levels.r] <- 0
        }
        if (any(g.levels.r) && is.element(struct[1], c("HCS", 
            "HAR", "DIAG"))) {
            G[g.levels.r, ] <- 0
            G[, g.levels.r] <- 0
            tau2[g.levels.r] <- 0
            warning("Fixed 'tau2' to 0 for removed level(s).")
        }
        if (any(g.levels.r) && struct[1] == "UN") {
            G[g.levels.r, ] <- 0
            G[, g.levels.r] <- 0
            tau2[g.levels.r] <- 0
            rho <- G[upper.tri(G)]
            warning("Fixed 'tau2' and corresponding 'rho' value(s) to 0 for removed level(s).")
        }
        if (any(g.levels.r) && struct[1] == "UNHO") {
            G[g.levels.r, ] <- 0
            G[, g.levels.r] <- 0
            diag(G) <- tau2
            rho <- G[upper.tri(G)]
            warning("Fixed 'rho' value(s) to 0 corresponding to removed level(s).")
        }
        if (g.nlevels.f[1] == 2) {
            if (is.element(struct[1], c("CS", "AR", "UNHO")) && 
                !is.na(tau2) && tau2 == 0) 
                rho <- 0
            if (is.element(struct[1], c("HCS", "UN", "HAR")) && 
                ((!is.na(tau2[1]) && tau2[1] == 0) || (!is.na(tau2[2]) && 
                  tau2[2] == 0))) 
                rho <- 0
        }
    }
    else {
        G <- NULL
        g.levels.f <- NULL
        g.levels.r <- NULL
        g.levels.k <- NULL
        g.levels.comb.k <- NULL
        g.nlevels.f <- NULL
    }
    if (withH) {
        h.nlevels.f <- h.nlevels
        h.levels.f <- h.levels
        mf.h <- data.frame(inner = factor(mf.h[[1]]), outer = factor(mf.h[[2]]))
        h.nlevels <- c(nlevels(mf.h[[1]]), nlevels(mf.h[[2]]))
        h.levels <- list(levels(mf.h[[1]]), levels(mf.h[[2]]))
        h.levels.r <- !is.element(h.levels.f[[1]], h.levels[[1]])
        if (any(h.levels.r)) 
            warning("One or more levels of inner factor removed due to NAs.")
        if (is.element(struct[2], c("ID", "DIAG"))) 
            phi <- 0
        if (h.nlevels[1] == 1 && is.element(struct[2], c("CS", 
            "HCS", "AR", "HAR")) && is.na(phi)) {
            phi <- 0
            warning("Inner factor has only a single level, so fixed value of 'phi' to 0.")
        }
        h.levels.k <- table(factor(mf.h[[1]], levels = h.levels.f[[1]]))
        if (is.element(struct[2], c("HCS", "UN", "DIAG", "HAR"))) {
            if (any(is.na(gamma2) & h.levels.k == 1)) {
                gamma2[is.na(gamma2) & h.levels.k == 1] <- 0
                warning("Inner factor has k=1 for one or more levels. Corresponding gamma2 value(s) fixed to 0.")
            }
        }
        h.levels.comb.k <- crossprod(Z.H2, Z.H1)
        h.levels.comb.k <- split(h.levels.comb.k, 1:nrow(h.levels.comb.k))
        if (all(unlist(lapply(h.levels.comb.k, sum)) == 1)) {
            if (is.element(struct[2], c("CS", "HCS", "AR", "HAR")) && 
                is.na(phi)) {
                phi <- 0
                warning("Each level of the outer factor contains only a single level of the inner factor, so fixed value of 'phi' to 0.")
            }
        }
        h.levels.comb.k <- lapply(h.levels.comb.k, function(x) outer(x, 
            x, FUN = "&"))
        h.levels.comb.k <- lapply(h.levels.comb.k, function(x) ifelse(x, 
            1, 0))
        h.levels.comb.k <- Reduce("+", h.levels.comb.k)
        h.levels.comb.k <- h.levels.comb.k[upper.tri(h.levels.comb.k)]
        if (is.element(struct[2], c("UN", "UNHO")) && any(h.levels.comb.k == 
            0 & is.na(phi))) {
            phi[h.levels.comb.k == 0] <- 0
            warning("Some combinations of the levels of the inner factor never occurred. Corresponding 'phi' value(s) fixed to 0.")
        }
        if (is.element(struct[2], c("UN", "UNHO")) && h.nlevels.f[1] == 
            1 && is.na(phi)) {
            phi <- 0
            warning("Inner factor has only a single level, so fixed value of 'phi' to 0.")
        }
        if (struct[2] == "CS") {
            H <- matrix(phi * gamma2, nrow = h.nlevels.f[1], 
                ncol = h.nlevels.f[1])
            diag(H) <- gamma2
        }
        if (struct[2] == "HCS") {
            H <- matrix(phi, nrow = h.nlevels.f[1], ncol = h.nlevels.f[1])
            diag(H) <- 1
            H <- diag(sqrt(gamma2), nrow = h.nlevels.f[1], ncol = h.nlevels.f[1]) %*% 
                H %*% diag(sqrt(gamma2), nrow = h.nlevels.f[1], 
                ncol = h.nlevels.f[1])
            diag(H) <- gamma2
        }
        if (struct[2] == "UN") {
            H <- .con.vcov.UN(gamma2, phi)
        }
        if (struct[2] == "ID" || struct[2] == "DIAG") {
            H <- diag(gamma2, nrow = h.nlevels.f[1], ncol = h.nlevels.f[1])
        }
        if (struct[2] == "UNHO") {
            H <- matrix(NA_real_, nrow = h.nlevels.f[1], ncol = h.nlevels.f[1])
            H[upper.tri(H)] <- phi
            H[lower.tri(H)] <- t(H)[lower.tri(H)]
            diag(H) <- 1
            H <- diag(sqrt(rep(gamma2, h.nlevels.f[1])), nrow = h.nlevels.f[1], 
                ncol = h.nlevels.f[1]) %*% H %*% diag(sqrt(rep(gamma2, 
                h.nlevels.f[1])), nrow = h.nlevels.f[1], ncol = h.nlevels.f[1])
            diag(H) <- gamma2
        }
        if (struct[2] == "AR") {
            if (is.na(phi)) {
                H <- matrix(NA_real_, nrow = h.nlevels.f[1], 
                  ncol = h.nlevels.f[1])
            }
            else {
                if (h.nlevels.f[1] > 1) {
                  H <- toeplitz(ARMAacf(ar = phi, lag.max = h.nlevels.f[1] - 
                    1))
                }
                else {
                  H <- diag(1)
                }
            }
            H <- diag(sqrt(rep(gamma2, h.nlevels.f[1])), nrow = h.nlevels.f[1], 
                ncol = h.nlevels.f[1]) %*% H %*% diag(sqrt(rep(gamma2, 
                h.nlevels.f[1])), nrow = h.nlevels.f[1], ncol = h.nlevels.f[1])
            diag(H) <- gamma2
        }
        if (struct[2] == "HAR") {
            if (is.na(phi)) {
                H <- matrix(NA_real_, nrow = h.nlevels.f[1], 
                  ncol = h.nlevels.f[1])
            }
            else {
                if (h.nlevels.f[1] > 1) {
                  H <- toeplitz(ARMAacf(ar = phi, lag.max = h.nlevels.f[1] - 
                    1))
                }
                else {
                  H <- diag(1)
                }
            }
            H <- diag(sqrt(gamma2), nrow = h.nlevels.f[1], ncol = h.nlevels.f[1]) %*% 
                H %*% diag(sqrt(gamma2), nrow = h.nlevels.f[1], 
                ncol = h.nlevels.f[1])
            diag(H) <- gamma2
        }
        if (any(h.levels.r) && is.element(struct[2], c("CS", 
            "AR", "ID"))) {
            H[h.levels.r, ] <- 0
            H[, h.levels.r] <- 0
        }
        if (any(h.levels.r) && is.element(struct[2], c("HCS", 
            "HAR", "DIAG"))) {
            H[h.levels.r, ] <- 0
            H[, h.levels.r] <- 0
            gamma2[h.levels.r] <- 0
            warning("Fixed 'gamma2' to 0 for removed level(s).")
        }
        if (any(h.levels.r) && struct[2] == "UN") {
            H[h.levels.r, ] <- 0
            H[, h.levels.r] <- 0
            gamma2[h.levels.r] <- 0
            phi <- H[upper.tri(H)]
            warning("Fixed 'gamma2' and corresponding 'phi' value(s) to 0 for removed level(s).")
        }
        if (any(h.levels.r) && struct[2] == "UNHO") {
            H[h.levels.r, ] <- 0
            H[, h.levels.r] <- 0
            diag(H) <- gamma2
            phi <- H[upper.tri(H)]
            warning("Fixed 'phi' value(s) to 0 corresponding to removed level(s).")
        }
        if (h.nlevels.f[1] == 2) {
            if (is.element(struct[2], c("CS", "AR", "UNHO")) && 
                !is.na(gamma2) && gamma2 == 0) 
                phi <- 0
            if (is.element(struct[2], c("HCS", "UN", "HAR")) && 
                ((!is.na(gamma2[1]) && gamma2[1] == 0) || (!is.na(gamma2[2]) && 
                  gamma2[2] == 0))) 
                phi <- 0
        }
    }
    else {
        H <- NULL
        h.levels.f <- NULL
        h.levels.r <- NULL
        h.levels.k <- NULL
        h.levels.comb.k <- NULL
        h.nlevels.f <- NULL
    }
    Y <- as.matrix(yi)
    if (very.verbose) 
        message("Extracting/computing initial values ...")
    if (verbose) {
        L.FE <- try(chol(V), silent = TRUE)
    }
    else {
        L.FE <- suppressWarnings(try(chol(V), silent = !verbose))
    }
    if (inherits(L.FE, "try-error")) {
        sigma2.init <- rep(0.001, sigma2s)
        tau2.init <- rep(0.001, tau2s)
        rho.init <- rep(0.5, rhos)
        gamma2.init <- rep(0.001, gamma2s)
        phi.init <- rep(0.5, rhos)
    }
    else {
        W <- chol2inv(L.FE)
        U <- chol(W)
        sX <- U %*% X
        sY <- U %*% Y
        b.FE <- solve(crossprod(sX), crossprod(sX, sY))
        total <- max(0.001 * (sigma2s + tau2s + gamma2s), var(as.vector(Y) - 
            as.vector(X %*% b.FE)) - 1/mean(1/diag(V)))
        sigma2.init <- rep(total/(sigma2s + tau2s + gamma2s), 
            sigma2s)
        tau2.init <- rep(total/(sigma2s + tau2s + gamma2s), tau2s)
        gamma2.init <- rep(total/(sigma2s + tau2s + gamma2s), 
            gamma2s)
        rho.init <- rep(0.5, rhos)
        phi.init <- rep(0.5, phis)
        QE <- sum(as.vector(sY - sX %*% b.FE)^2)
    }
    con <- list(verbose = FALSE, optimizer = "nlminb", optmethod = "BFGS", 
        sigma2.init = sigma2.init, tau2.init = tau2.init, rho.init = rho.init, 
        gamma2.init = gamma2.init, phi.init = phi.init, REMLf = TRUE, 
        tol = 1e-07, cholesky = ifelse(struct == "UN", TRUE, 
            FALSE), posdefify = FALSE, hessian = FALSE, vctransf = FALSE)
    con.pos <- pmatch(names(control), names(con))
    con[c(na.omit(con.pos))] <- control[!is.na(con.pos)]
    if (verbose) 
        con$verbose <- verbose
    verbose <- con$verbose
    if (withS && any(con$sigma2.init <= 0)) 
        stop("Values of 'sigma2.init' must be positive.")
    if (withG && any(con$tau2.init <= 0)) 
        stop("Values of 'tau2.init' must be positive.")
    if (withG && any(con$rho.init <= -1 | con$rho.init >= 1)) 
        stop("Values of 'rho.init' must be in (-1,1).")
    if (withH && any(con$gamma2.init <= 0)) 
        stop("Values of 'gamma2.init' must be positive.")
    if (withH && any(con$phi.init <= -1 | con$phi.init >= 1)) 
        stop("Values of 'phi.init' must be in (-1,1).")
    if (!withG) 
        con$cholesky[1] <- FALSE
    if (con$cholesky[1] && struct[1] != "UN") 
        con$cholesky[1] <- FALSE
    if (!withH) 
        con$cholesky[2] <- FALSE
    if (con$cholesky[2] && struct[2] != "UN") 
        con$cholesky[2] <- FALSE
    sigma2.init <- con$sigma2.init
    tau2.init <- con$tau2.init
    rho.init <- con$rho.init
    gamma2.init <- con$gamma2.init
    phi.init <- con$phi.init
    con$sigma2.init <- log(sigma2.init)
    if (con$cholesky[1]) {
        G <- .con.vcov.UN(tau2.init, rho.init)
        G <- try(chol(G), silent = TRUE)
        if (inherits(G, "try-error")) 
            stop("Cannot take Choleski decomposition of initial G matrix.")
        con$tau2.init <- diag(G)
        con$rho.init <- G[upper.tri(G)]
    }
    else {
        con$tau2.init <- log(tau2.init)
        con$rho.init <- transf.rtoz(rho.init)
    }
    if (con$cholesky[2]) {
        H <- .con.vcov.UN(gamma2.init, phi.init)
        H <- try(chol(H), silent = TRUE)
        if (inherits(H, "try-error")) 
            stop("Cannot take Choleski decomposition of initial H matrix.")
        con$gamma2.init <- diag(H)
        con$phi.init <- H[upper.tri(H)]
    }
    else {
        con$gamma2.init <- log(gamma2.init)
        con$phi.init <- transf.rtoz(phi.init)
    }
    optimizer <- match.arg(con$optimizer, c("optim", "nlminb", 
        "uobyqa", "newuoa", "bobyqa", "nloptr", "nlm", "hjk", 
        "nmk", "ucminf"))
    optmethod <- con$optmethod
    tol <- con$tol
    posdefify <- con$posdefify
    cholesky <- con$cholesky
    optcontrol <- control[is.na(con.pos)]
    if (length(optcontrol) == 0) 
        optcontrol <- list()
    if (optimizer == "nloptr" && !is.element("algorithm", names(optcontrol))) 
        optcontrol$algorithm <- "NLOPT_LN_BOBYQA"
    if (optimizer == "nloptr" && !is.element("ftol_rel", names(optcontrol))) 
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
    if (is.element(optimizer, c("hjk", "nmk"))) {
        if (!requireNamespace("dfoptim", quietly = TRUE)) 
            stop("Please install the 'dfoptim' package to use this optimizer.")
    }
    if (optimizer == "ucminf") {
        if (!requireNamespace("ucminf", quietly = TRUE)) 
            stop("Please install the 'ucminf' package to use this optimizer.")
    }
    if (withS) {
        if (length(con$sigma2.init) != sigma2s) 
            stop(paste("Length of 'sigma2.init' argument (", 
                length(con$sigma2.init), ") does not match actual number of variance components (", 
                sigma2s, ").", sep = ""))
    }
    else {
        con$sigma2.init <- 0
    }
    if (withG) {
        if (length(con$tau2.init) != tau2s) 
            stop(paste("Length of 'tau2.init' argument (", length(con$tau2.init), 
                ") does not match actual number of variance components (", 
                tau2s, ").", sep = ""))
    }
    else {
        con$tau2.init <- 0
    }
    if (withG) {
        if (length(con$rho.init) != rhos) 
            stop(paste("Length of 'rho.init' argument (", length(con$rho.init), 
                ") does not match actual number of correlations (", 
                rhos, ").", sep = ""))
    }
    else {
        con$rho.init <- 0
    }
    if (withH) {
        if (length(con$gamma2.init) != gamma2s) 
            stop(paste("Length of 'gamma2.init' argument (", 
                length(con$gamma2.init), ") does not match actual number of variance components (", 
                gamma2s, ").", sep = ""))
    }
    else {
        con$gamma2.init <- 0
    }
    if (withH) {
        if (length(con$phi.init) != phis) 
            stop(paste("Length of 'phi.init' argument (", length(con$phi.init), 
                ") does not match actual number of correlations (", 
                phis, ").", sep = ""))
    }
    else {
        con$phi.init <- 0
    }
    if (any(eigen(crossprod(X), symmetric = TRUE, only.values = TRUE)$values <= 
        tol)) 
        stop("Model matrix not of full rank. Cannot fit model.")
    if (withS) {
        sigma2.fix <- !is.na(sigma2)
    }
    else {
        sigma2.fix <- NA
    }
    if (withG) {
        tau2.fix <- !is.na(tau2)
        rho.fix <- !is.na(rho)
    }
    else {
        tau2.fix <- NA
        rho.fix <- NA
    }
    if (withH) {
        gamma2.fix <- !is.na(gamma2)
        phi.fix <- !is.na(phi)
    }
    else {
        gamma2.fix <- NA
        phi.fix <- NA
    }
    vc.fix <- list(sigma2 = sigma2.fix, tau2 = tau2.fix, rho = rho.fix, 
        gamma2 = gamma2.fix, phi = phi.fix)
    if (verbose) {
        cat("\nVariance Components in Model:")
        if (!withS && !withG && !withH) {
            cat(" none\n\n")
        }
        else {
            cat("\n\n")
            vcs <- rbind(c(sigma2 = if (withS) round(sigma2.init, 
                digits = digits) else NA, tau2 = if (withG) round(tau2.init, 
                digits = digits) else NA, rho = if (withG) round(rho.init, 
                digits = digits) else NA, gamma2 = if (withH) round(gamma2.init, 
                digits = digits) else NA, phi = if (withH) round(phi.init, 
                digits = digits) else NA), round(c(if (withS) sigma2 else NA, 
                if (withG) tau2 else NA, if (withG) rho else NA, 
                if (withH) gamma2 else NA, if (withH) phi else NA), 
                digits = digits))
            vcs <- data.frame(vcs)
            rownames(vcs) <- c("initial", "specified")
            vcs <- rbind(included = ifelse(c(rep(withS, sigma2s), 
                rep(withG, tau2s), rep(withG, rhos), rep(withH, 
                  gamma2s), rep(withH, phis)), "Yes", "No"), 
                fixed = unlist(vc.fix), vcs)
            print(vcs, na.print = "")
            cat("\n")
        }
    }
    alpha <- ifelse(level > 1, (100 - level)/100, 1 - level)
    if (very.verbose) 
        message("Model fitting ...")
    if (optimizer == "optim") {
        par.arg <- "par"
        ctrl.arg <- ", control=optcontrol"
    }
    if (optimizer == "nlminb") {
        par.arg <- "start"
        ctrl.arg <- ", control=optcontrol"
    }
    if (is.element(optimizer, c("uobyqa", "newuoa", "bobyqa"))) {
        par.arg <- "par"
        optimizer <- paste0("minqa::", optimizer)
        ctrl.arg <- ", control=optcontrol"
    }
    if (optimizer == "nloptr") {
        par.arg <- "x0"
        optimizer <- paste0("nloptr::nloptr")
        ctrl.arg <- ", opts=optcontrol"
    }
    if (optimizer == "nlm") {
        par.arg <- "p"
        ctrl.arg <- paste(names(optcontrol), unlist(optcontrol), 
            sep = "=", collapse = ", ")
        if (nchar(ctrl.arg) != 0) 
            ctrl.arg <- paste0(", ", ctrl.arg)
    }
    if (is.element(optimizer, c("hjk", "nmk"))) {
        par.arg <- "par"
        optimizer <- paste0("dfoptim::", optimizer)
        ctrl.arg <- ", control=optcontrol"
    }
    if (optimizer == "ucminf") {
        par.arg <- "par"
        optimizer <- paste0("ucminf::ucminf")
        ctrl.arg <- ", control=optcontrol"
    }
    optcall <- paste(optimizer, "(", par.arg, "=c(con$sigma2.init, con$tau2.init, con$rho.init, con$gamma2.init, con$phi.init),\n      .ll.rma.mv, reml=reml, ", 
        ifelse(optimizer == "optim", "method=optmethod, ", ""), 
        "Y=Y, M=V, X.fit=X, k=k, pX=p,\n      D.S=D.S, Z.G1=Z.G1, Z.G2=Z.G2, Z.H1=Z.H1, Z.H2=Z.H2,\n      sigma2.val=sigma2, tau2.val=tau2, rho.val=rho, gamma2.val=gamma2, phi.val=phi,\n      sigma2s=sigma2s, tau2s=tau2s, rhos=rhos, gamma2s=gamma2s, phis=phis,\n      withS=withS, withG=withG, withH=withH,\n      struct=struct, g.levels.r=g.levels.r, h.levels.r=h.levels.r,\n      tol=tol, sparse=sparse, cholesky=cholesky, posdefify=posdefify, vctransf=TRUE,\n      verbose=verbose, very.verbose=very.verbose, digits=digits, REMLf=con$REMLf", 
        ctrl.arg, ")\n", sep = "")
    opt.res <- try(eval(parse(text = optcall)), silent = !verbose)
    if (inherits(opt.res, "try-error")) 
        stop("Error during optimization.")
    if (is.element(optimizer, c("optim", "nlminb", "dfoptim::hjk", 
        "dfoptim::nmk")) && opt.res$convergence != 0) 
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
    if (optimizer == "ucminf::ucminf" && !(opt.res$convergence == 
        1 || opt.res$convergence == 2)) 
        stop(paste0("Optimizer (", optimizer, ") did not achieve convergence (convergence = ", 
            opt.res$convergence, ")."))
    if (is.element(optimizer, c("optim", "dfoptim::hjk", "dfoptim::nmk", 
        "ucminf::ucminf"))) 
        ll <- -1 * c(opt.res$value)
    if (is.element(optimizer, c("nlminb", "nloptr::nloptr"))) 
        ll <- -1 * c(opt.res$objective)
    if (is.element(optimizer, c("minqa::uobyqa", "minqa::newuoa", 
        "minqa::bobyqa"))) 
        ll <- -1 * c(opt.res$fval)
    if (optimizer == "nlm") 
        ll <- -1 * c(opt.res$minimum)
    if (verbose) {
        cat("\n")
        print(opt.res)
    }
    if (optimizer == "nloptr::nloptr") 
        opt.res$par <- opt.res$solution
    if (optimizer == "nlm") 
        opt.res$par <- opt.res$estimate
    if (p < k) {
        sigma2.val <- sigma2
        tau2.val <- tau2
        rho.val <- rho
        gamma2.val <- gamma2
        phi.val <- phi
        sigma2 <- ifelse(is.na(sigma2), exp(opt.res$par[1:sigma2s]), 
            sigma2)
        if (cholesky[1]) {
            tau2 <- opt.res$par[(sigma2s + 1):(sigma2s + tau2s)]
            rho <- opt.res$par[(sigma2s + tau2s + 1):(sigma2s + 
                tau2s + rhos)]
        }
        else {
            tau2 <- ifelse(is.na(tau2), exp(opt.res$par[(sigma2s + 
                1):(sigma2s + tau2s)]), tau2)
            rho <- ifelse(is.na(rho), transf.ztor(opt.res$par[(sigma2s + 
                tau2s + 1):(sigma2s + tau2s + rhos)]), rho)
            tau2 <- ifelse(tau2 <= .Machine$double.eps * 10, 
                0, tau2)
        }
        if (cholesky[2]) {
            gamma2 <- opt.res$par[(sigma2s + tau2s + rhos + 1):(sigma2s + 
                tau2s + rhos + gamma2s)]
            phi <- opt.res$par[(sigma2s + tau2s + rhos + gamma2s + 
                1):(sigma2s + tau2s + rhos + gamma2s + phis)]
        }
        else {
            gamma2 <- ifelse(is.na(gamma2), exp(opt.res$par[(sigma2s + 
                tau2s + rhos + 1):(sigma2s + tau2s + rhos + gamma2s)]), 
                gamma2)
            phi <- ifelse(is.na(phi), transf.ztor(opt.res$par[(sigma2s + 
                tau2s + rhos + gamma2s + 1):(sigma2s + tau2s + 
                rhos + gamma2s + phis)]), phi)
            gamma2 <- ifelse(gamma2 <= .Machine$double.eps * 
                10, 0, gamma2)
        }
        sigma2 <- ifelse(sigma2 <= .Machine$double.eps * 10, 
            0, sigma2)
    }
    else {
        sigma2[is.na(sigma2)] <- 0
        tau2[is.na(tau2)] <- 0
        rho[is.na(rho)] <- 0
        gamma2[is.na(gamma2)] <- 0
        phi[is.na(phi)] <- 0
    }
    M <- V
    if (withS) {
        for (j in seq_len(sigma2s)) {
            M <- M + sigma2[j] * D.S[[j]]
        }
    }
    if (withG) {
        ncol.Z.G1 <- ncol(Z.G1)
        if (struct[1] == "CS") {
            G <- matrix(rho * tau2, nrow = ncol.Z.G1, ncol = ncol.Z.G1)
            diag(G) <- tau2
        }
        if (struct[1] == "HCS") {
            G <- matrix(rho, nrow = ncol.Z.G1, ncol = ncol.Z.G1)
            diag(G) <- 1
            G <- diag(sqrt(tau2), nrow = ncol.Z.G1, ncol = ncol.Z.G1) %*% 
                G %*% diag(sqrt(tau2), nrow = ncol.Z.G1, ncol = ncol.Z.G1)
            diag(G) <- tau2
        }
        if (struct[1] == "UN") {
            if (cholesky[1]) {
                G <- .con.vcov.UN.chol(tau2, rho)
                tau2 <- diag(G)
                rho <- cov2cor(G)[upper.tri(G)]
                tau2[!is.na(tau2.val)] <- tau2.val[!is.na(tau2.val)]
                rho[!is.na(rho.val)] <- rho.val[!is.na(rho.val)]
            }
            G <- .con.vcov.UN(tau2, rho)
            if (posdefify) {
                G <- as.matrix(nearPD(G)$mat)
                tau2 <- diag(G)
                rho <- cov2cor(G)[upper.tri(G)]
            }
        }
        if (struct[1] == "ID" || struct[1] == "DIAG") {
            G <- diag(tau2, nrow = ncol.Z.G1, ncol = ncol.Z.G1)
        }
        if (struct[1] == "UNHO") {
            G <- matrix(NA_real_, nrow = ncol.Z.G1, ncol = ncol.Z.G1)
            G[upper.tri(G)] <- rho
            G[lower.tri(G)] <- t(G)[lower.tri(G)]
            diag(G) <- 1
            G <- diag(sqrt(rep(tau2, ncol.Z.G1)), nrow = ncol.Z.G1, 
                ncol = ncol.Z.G1) %*% G %*% diag(sqrt(rep(tau2, 
                ncol.Z.G1)), nrow = ncol.Z.G1, ncol = ncol.Z.G1)
            if (posdefify) {
                G <- as.matrix(nearPD(G, keepDiag = TRUE)$mat)
                tau2 <- G[1, 1]
                rho <- cov2cor(G)[upper.tri(G)]
            }
        }
        if (struct[1] == "AR") {
            if (ncol.Z.G1 > 1) {
                G <- toeplitz(ARMAacf(ar = rho, lag.max = ncol.Z.G1 - 
                  1))
            }
            else {
                G <- diag(1)
            }
            G <- diag(sqrt(rep(tau2, ncol.Z.G1)), nrow = ncol.Z.G1, 
                ncol = ncol.Z.G1) %*% G %*% diag(sqrt(rep(tau2, 
                ncol.Z.G1)), nrow = ncol.Z.G1, ncol = ncol.Z.G1)
            diag(G) <- tau2
        }
        if (struct[1] == "HAR") {
            if (ncol.Z.G1 > 1) {
                G <- toeplitz(ARMAacf(ar = rho, lag.max = ncol.Z.G1 - 
                  1))
            }
            else {
                G <- diag(1)
            }
            G <- diag(sqrt(tau2), nrow = ncol.Z.G1, ncol = ncol.Z.G1) %*% 
                G %*% diag(sqrt(tau2), nrow = ncol.Z.G1, ncol = ncol.Z.G1)
            diag(G) <- tau2
        }
        if (any(g.levels.r)) {
            G[g.levels.r, ] <- 0
            G[, g.levels.r] <- 0
        }
        if (sparse) 
            G <- Matrix(G, sparse = TRUE)
        M <- M + (Z.G1 %*% G %*% t(Z.G1)) * tcrossprod(Z.G2)
        colnames(G) <- rownames(G) <- g.levels.f[[1]]
    }
    if (withH) {
        ncol.Z.H1 <- ncol(Z.H1)
        if (struct[2] == "CS") {
            H <- matrix(phi * gamma2, nrow = ncol.Z.H1, ncol = ncol.Z.H1)
            diag(H) <- gamma2
        }
        if (struct[2] == "HCS") {
            H <- matrix(phi, nrow = ncol.Z.H1, ncol = ncol.Z.H1)
            diag(H) <- 1
            H <- diag(sqrt(gamma2), nrow = ncol.Z.H1, ncol = ncol.Z.H1) %*% 
                H %*% diag(sqrt(gamma2), nrow = ncol.Z.H1, ncol = ncol.Z.H1)
            diag(H) <- gamma2
        }
        if (struct[2] == "UN") {
            if (cholesky[2]) {
                H <- .con.vcov.UN.chol(gamma2, phi)
                gamma2 <- diag(H)
                phi <- cov2cor(H)[upper.tri(H)]
                gamma2[!is.na(gamma2.val)] <- gamma2.val[!is.na(gamma2.val)]
                phi[!is.na(phi.val)] <- phi.val[!is.na(phi.val)]
            }
            H <- .con.vcov.UN(gamma2, phi)
            if (posdefify) {
                H <- as.matrix(nearPD(H)$mat)
                gamma2 <- diag(H)
                phi <- cov2cor(H)[upper.tri(H)]
            }
        }
        if (struct[2] == "ID" || struct[2] == "DIAG") {
            H <- diag(gamma2, nrow = ncol.Z.H1, ncol = ncol.Z.H1)
        }
        if (struct[2] == "UNHO") {
            H <- matrix(NA_real_, nrow = ncol.Z.H1, ncol = ncol.Z.H1)
            H[upper.tri(H)] <- phi
            H[lower.tri(H)] <- t(H)[lower.tri(H)]
            diag(H) <- 1
            H <- diag(sqrt(rep(gamma2, ncol.Z.H1)), nrow = ncol.Z.H1, 
                ncol = ncol.Z.H1) %*% H %*% diag(sqrt(rep(gamma2, 
                ncol.Z.H1)), nrow = ncol.Z.H1, ncol = ncol.Z.H1)
            if (posdefify) {
                H <- as.matrix(nearPD(H, keepDiag = TRUE)$mat)
                gamma2 <- H[1, 1]
                phi <- cov2cor(H)[upper.tri(H)]
            }
        }
        if (struct[2] == "AR") {
            if (ncol.Z.H1 > 1) {
                H <- toeplitz(ARMAacf(ar = phi, lag.max = ncol.Z.H1 - 
                  1))
            }
            else {
                H <- diag(1)
            }
            H <- diag(sqrt(rep(gamma2, ncol.Z.H1)), nrow = ncol.Z.H1, 
                ncol = ncol.Z.H1) %*% H %*% diag(sqrt(rep(gamma2, 
                ncol.Z.H1)), nrow = ncol.Z.H1, ncol = ncol.Z.H1)
            diag(H) <- gamma2
        }
        if (struct[2] == "HAR") {
            if (ncol.Z.H1 > 1) {
                H <- toeplitz(ARMAacf(ar = phi, lag.max = ncol.Z.H1 - 
                  1))
            }
            else {
                H <- diag(1)
            }
            H <- diag(sqrt(gamma2), nrow = ncol.Z.H1, ncol = ncol.Z.H1) %*% 
                H %*% diag(sqrt(gamma2), nrow = ncol.Z.H1, ncol = ncol.Z.H1)
            diag(H) <- gamma2
        }
        if (any(h.levels.r)) {
            H[h.levels.r, ] <- 0
            H[, h.levels.r] <- 0
        }
        if (sparse) 
            H <- Matrix(H, sparse = TRUE)
        M <- M + (Z.H1 %*% H %*% t(Z.H1)) * tcrossprod(Z.H2)
        colnames(H) <- rownames(H) <- h.levels.f[[1]]
    }
    colnames(M) <- rownames(M) <- NULL
    if (posdefify) 
        M <- as.matrix(nearPD(M)$mat)
    if (verbose) {
        L <- try(chol(M), silent = !verbose)
    }
    else {
        L <- suppressWarnings(try(chol(M), silent = !verbose))
    }
    if (inherits(L, "try-error")) {
        stop("Final variance-covariance matrix not positive definite.")
    }
    else {
        W <- chol2inv(L)
        U <- chol(W)
        sX <- U %*% X
    }
    if (is.null(A)) {
        sY <- U %*% Y
        vb <- matrix(solve(crossprod(sX)), nrow = p, ncol = p)
        b <- matrix(vb %*% crossprod(sX, sY), ncol = 1)
        RSS.f <- sum(as.vector(sY - sX %*% b)^2)
    }
    else {
        stXAX <- chol2inv(chol(as.matrix(t(X) %*% A %*% X)))
        b <- matrix(stXAX %*% crossprod(X, A) %*% Y, ncol = 1)
        vb <- matrix(stXAX %*% t(X) %*% A %*% M %*% A %*% X %*% 
            stXAX, nrow = p, ncol = p)
        RSS.f <- as.vector(t(Y - X %*% b) %*% W %*% (Y - X %*% 
            b))
    }
    rownames(vb) <- colnames(vb) <- rownames(b) <- colnames(X)
    se <- sqrt(diag(vb))
    names(se) <- NULL
    zval <- c(b/se)
    pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
    QM <- as.vector(t(b)[btt] %*% chol2inv(chol(vb[btt, btt])) %*% 
        b[btt])
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
    QE.df <- k - p
    if (QE.df > 0L) {
        if (inherits(L.FE, "try-error")) {
            QE <- NA
            QEp <- NA
        }
        else {
            QEp <- pchisq(QE, df = QE.df, lower.tail = FALSE)
        }
    }
    else {
        QE <- 0
        QEp <- 1
    }
    ll.QE <- -1/2 * (k) * log(2 * base::pi) - 1/2 * determinant(V, 
        logarithm = TRUE)$modulus
    if (con$hessian) {
        if (very.verbose) 
            message("Computing Hessian ...")
        if (!requireNamespace("numDeriv", quietly = TRUE)) 
            stop("Please install the 'numDeriv' package for Hessian computation.")
        hessian <- try(numDeriv::hessian(func = .ll.rma.mv, x = if (con$vctransf) 
            opt.res$par
        else c(sigma2, tau2, rho, gamma2, phi), reml = reml, 
            Y = Y, M = V, X.fit = X, k = k, pX = p, D.S = D.S, 
            Z.G1 = Z.G1, Z.G2 = Z.G2, Z.H1 = Z.H1, Z.H2 = Z.H2, 
            sigma2.val = sigma2.val, tau2.val = tau2.val, rho.val = rho.val, 
            gamma2.val = gamma2.val, phi.val = phi.val, sigma2s = sigma2s, 
            tau2s = tau2s, rhos = rhos, gamma2s = gamma2s, phis = phis, 
            withS = withS, withG = withG, withH = withH, struct = struct, 
            g.levels.r = g.levels.r, h.levels.r = h.levels.r, 
            tol = tol, sparse = sparse, cholesky = ifelse(c(con$vctransf, 
                con$vctransf) & cholesky, TRUE, FALSE), posdefify = posdefify, 
            vctransf = con$vctransf, verbose = verbose, very.verbose = very.verbose, 
            digits = digits, REMLf = con$REMLf), silent = TRUE)
        if (inherits(hessian, "try-error")) 
            warning("Error when trying to compute Hessian.")
        colnames(hessian) <- 1:ncol(hessian)
        if (length(sigma2) == 1) {
            colnames(hessian)[1] <- "sigma^2"
        }
        else {
            colnames(hessian)[1:length(sigma2)] <- paste("sigma^2.", 
                1:length(sigma2), sep = "")
        }
        if (length(tau2) == 1) {
            colnames(hessian)[length(sigma2) + 1] <- "tau^2"
        }
        else {
            colnames(hessian)[(length(sigma2) + 1):(length(sigma2) + 
                length(tau2))] <- paste("tau^2.", 1:length(tau2), 
                sep = "")
        }
        if (length(rho) == 1) {
            colnames(hessian)[length(sigma2) + length(tau2) + 
                1] <- "rho"
        }
        else {
            colnames(hessian)[(length(sigma2) + length(tau2) + 
                1):(length(sigma2) + length(tau2) + length(rho))] <- paste("rho.", 
                outer(1:g.nlevels.f[1], 1:g.nlevels.f, paste, 
                  sep = ".")[upper.tri(matrix(NA, nrow = g.nlevels.f, 
                  ncol = g.nlevels.f))], sep = "")
        }
        if (length(gamma2) == 1) {
            colnames(hessian)[length(sigma2) + length(tau2) + 
                length(rho) + 1] <- "gamma^2"
        }
        else {
            colnames(hessian)[(length(sigma2) + length(tau2) + 
                length(rho) + 1):(length(sigma2) + length(tau2) + 
                length(rho) + length(gamma2))] <- paste("gamma^2.", 
                1:length(gamma2), sep = "")
        }
        if (length(phi) == 1) {
            colnames(hessian)[length(sigma2) + length(tau2) + 
                length(rho) + length(gamma2) + 1] <- "phi"
        }
        else {
            colnames(hessian)[(length(sigma2) + length(tau2) + 
                length(rho) + length(gamma2) + 1):(length(sigma2) + 
                length(tau2) + length(rho) + length(gamma2) + 
                length(phi))] <- paste("phi.", outer(1:h.nlevels.f[1], 
                1:h.nlevels.f, paste, sep = ".")[upper.tri(matrix(NA, 
                nrow = h.nlevels.f, ncol = h.nlevels.f))], sep = "")
        }
        rownames(hessian) <- colnames(hessian)
        if (withS && withG && !withH) 
            hessian <- hessian[1:(nrow(hessian) - 2), 1:(ncol(hessian) - 
                2), drop = FALSE]
        if (withS && !withG && !withH) 
            hessian <- hessian[1:(nrow(hessian) - 4), 1:(ncol(hessian) - 
                4), drop = FALSE]
        if (!withS && withG && withH) 
            hessian <- hessian[2:nrow(hessian), 2:ncol(hessian), 
                drop = FALSE]
        if (!withS && withG && !withH) 
            hessian <- hessian[2:(nrow(hessian) - 2), 2:(ncol(hessian) - 
                2), drop = FALSE]
        if (!withS && !withG && !withH) 
            hessian <- NA
    }
    else {
        hessian <- NA
    }
    if (very.verbose) 
        message("Computing fit statistics and log likelihood ...")
    parms <- p + ifelse(withS, sum(ifelse(sigma2.fix, 0, 1)), 
        0) + ifelse(withG, sum(ifelse(tau2.fix, 0, 1)), 0) + 
        ifelse(withG, sum(ifelse(rho.fix, 0, 1)), 0) + ifelse(withH, 
        sum(ifelse(gamma2.fix, 0, 1)), 0) + ifelse(withH, sum(ifelse(phi.fix, 
        0, 1)), 0)
    ll.ML <- -1/2 * (k) * log(2 * base::pi) - 1/2 * determinant(M, 
        logarithm = TRUE)$modulus - 1/2 * RSS.f
    ll.REML <- -1/2 * (k - p) * log(2 * base::pi) + ifelse(con$REMLf, 
        1/2 * determinant(crossprod(X), logarithm = TRUE)$modulus, 
        0) - 1/2 * determinant(M, logarithm = TRUE)$modulus - 
        1/2 * determinant(crossprod(X, W) %*% X, logarithm = TRUE)$modulus - 
        1/2 * RSS.f
    dev.ML <- -2 * (ll.ML - ll.QE)
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
    weighted <- TRUE
    res <- list(b = b, se = se, zval = zval, pval = pval, ci.lb = ci.lb, 
        ci.ub = ci.ub, vb = vb, sigma2 = sigma2, tau2 = tau2, 
        rho = rho, gamma2 = gamma2, phi = phi, k = k, k.f = k.f, 
        k.eff = k.eff, p = p, p.eff = p.eff, parms = parms, m = m, 
        QE = QE, QEp = QEp, QM = QM, QMp = QMp, int.only = int.only, 
        int.incl = int.incl, allvipos = allvipos, yi = yi, vi = vi, 
        V = V, W = A, X = X, yi.f = yi.f, vi.f = vi.f, V.f = V.f, 
        X.f = X.f, ni = ni, ni.f = ni.f, M = M, G = G, H = H, 
        hessian = hessian, ids = ids, not.na = not.na, subset = subset, 
        slab = slab, slab.null = slab.null, measure = measure, 
        method = method, weighted = weighted, knha = knha, dfs = dfs, 
        btt = btt, intercept = intercept, digits = digits, level = level, 
        sparse = sparse, control = control, fit.stats = fit.stats, 
        vc.fix = vc.fix, withS = withS, withG = withG, withH = withH, 
        withR = withR, sigma2s = sigma2s, tau2s = tau2s, rhos = rhos, 
        gamma2s = gamma2s, phis = phis, s.names = s.names, g.names = g.names, 
        h.names = h.names, s.nlevels = s.nlevels, g.nlevels.f = g.nlevels.f, 
        g.nlevels = g.nlevels, h.nlevels.f = h.nlevels.f, h.nlevels = h.nlevels, 
        g.levels.f = g.levels.f, g.levels.k = g.levels.k, g.levels.comb.k = g.levels.comb.k, 
        h.levels.f = h.levels.f, h.levels.k = h.levels.k, h.levels.comb.k = h.levels.comb.k, 
        struct = struct, Rfix = Rfix, R = R, Rscale = Rscale, 
        mf.r = mf.r, mf.g.f = mf.g.f, mf.h.f = mf.h.f, random = random, 
        version = packageVersion("metafor"), call = mf)
    class(res) <- c("rma.mv", "rma")
    return(res)
}
