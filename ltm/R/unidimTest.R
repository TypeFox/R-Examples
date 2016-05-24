unidimTest <-
function (object, data, thetas, IRT = TRUE, z.vals = NULL, B = 100, ...) {
    if (missing(object) && (missing(data) & missing(thetas)))
        stop("either 'object' or both 'data' and 'thetas' must be supplied.\n")
    if (missing(object)) {
        if (!inherits(data, "matrix") && !inherits(data, "data.frame"))
            stop("'data' must be either a data.frame or a matrix")
        data <- data.matrix(data)
        if (any(its <- apply(data, 2, function (x) { x <- x[!is.na(x)]; length(unique(x)) } ) > 2))
            stop("'data' contain more that 2 distinct values for item(s): ", paste(which(its), collapse = ", "))
        data <- apply(data, 2, function (x) if (all(unique(x) %in% c(1, 0, NA))) x else x - 1)
        data <- data[complete.cases(data), ]
        n <- nrow(data)
        if (n == 0)
            stop("zero rows after case-wise missing values deletion.\n")
        p <- ncol(data)
        if (nrow(thetas) != p)
            stop("the dimensions of 'data' and 'thetas' do not much.\n")
        parms <- thetas
    } else {
        if (!class(object) %in% c("ltm", "rasch", "tpm"))
            stop("Use only with 'ltm', 'rasch' or 'tpm' objects.\n")
        if (inherits(object, "ltm") && any(c(object$ltst$factors > 1, object$ltst$quad.z1)))
            stop("\nfor 'ltm' objects it is assumed that the two-parameter logistic model has been fitted\n\t(i.e., one latent variable and no nonlinear terms).")
        data <- object$X    
        data <- data.matrix(data)
        data <- data[complete.cases(data), ]
        n <- nrow(data)
        if (n == 0)
            stop("\nzero rows after case-wise missing values deletion.\n")
        p <- ncol(data)
        parms <- if (inherits(object, "tpm")) cbind(object$coef[, 2:3], plogis(object$coef[, 1])) else object$coef
        fsc <- factor.scores(object, resp.patterns = data)$score.dat
        ablts <- fsc$z1
        se.ablts <- fsc$se.z1
        IRT  <- FALSE
    }
    ind <- t(combn(p, 2))
    n.ind <- nrow(ind)
    eigenRho <- function (data, ...) {
        rho <- diag(p)
        for (i in 1:n.ind) {
            r <- ind[i, ]
            rho[rbind(r, rev(r))] <- polychor(data[, r[1]], data[, r[2]], ...)
        }
        rho. <- rho
        diag(rho.) <- rep(max(rho[ind]), p)
        list(Rho = rho, ev = eigen(rho., symmetric = TRUE, only.values = TRUE)$values)
    }
    eigR <- eigenRho(data)
    rho <- eigR$Rho
    Tobs <- eigR$ev
    T.boot <- matrix(0, B, length(Tobs))
    for (b in 1:B) {
        if (!missing(object))
            z.vals <- rnorm(n, ablts, se.ablts)
        data.new <- rmvlogis(n, parms, IRT = IRT, z.vals = z.vals)
        T.boot[b, ] <- eigenRho(data.new)$ev
    }
    pval <- (sum(T.boot[, 2] >= Tobs[2], na.rm = TRUE) + 1) / (B + 1)
    if (!is.null(cnams <- colnames(data)))
        dimnames(rho) <- list(cnams, cnams)
    out <- list(Tobs = Tobs, T.boot = T.boot, p.value = pval, Rho = rho, 
                call = if (missing(object)) NULL else object$call)
    class(out) <- "unidimTest"
    out
}
