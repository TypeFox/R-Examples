rdaFit <- function(X, Y, Z, scale = FALSE, ...) {
    sdFun <- function(m) {
        sqrt(colSums(sweep(m, 2, colMeans(m), "-")^2) / (nrow(m) - 1))
    }
    ZERO <- 1e-04
    X <- data.matrix(X)
    NR <- nrow(X) - 1
    Xbar <- scale(X, center = TRUE, scale = scale)
    SD <- sdFun(Xbar)
    if(scale)
        Xbar[is.nan(Xbar)] <- 0
    ##tot.chi <- sum(svd(Xbar, nu = 0, nv = 0)$d^2) / NR
    ## partial RDA?
    if(!missing(Z) && !is.null(Z)) {
        Z <- data.matrix(Z)
        Zr <- scale(Z, center = TRUE, scale = FALSE)
        Q <- qr(Zr)
        Z <- qr.fitted(Q, Xbar)
        if(zrank <- Q$rank)
            Xbar <- qr.resid(Q, Xbar)
    } else {
        Zr <- NULL
    }
    ## Do RDA
    Y <- data.matrix(Y)
    Yr <- scale(Y, center = TRUE, scale = FALSE)
    Q <- qr(cbind(Zr, Yr), tol = ZERO)
    if(is.null(Zr))
        yrank <- Q$rank
    else
        yrank <- Q$rank - zrank
    Y <- qr.fitted(Q, Xbar)
    sol <- svd(Y)
    rank <- min(yrank, sum(sol$d > ZERO))
    sol$d <- sol$d / sqrt(NR)
    ax.names <- paste("RDA", 1:length(sol$d), sep = "")
    colnames(sol$u) <- ax.names
    colnames(sol$v) <- ax.names
    names(sol$d) <- ax.names
    rownames(sol$u) <- rownames(X)
    rownames(sol$v) <- colnames(X)
    if(rank) {
        wa.eig <- (Xbar %*% sol$v[, 1:rank, drop = FALSE]) / sqrt(NR)
        RDA <- list(lambda = sol$d[1:rank]^2,
                    u = as.matrix(sol$u)[, 1:rank, drop = FALSE],
                    v = as.matrix(sol$v)[, 1:rank, drop = FALSE],
                    wa = sweep(wa.eig, 2, 1/sol$d[1:rank], "*"))
        oo <- Q$pivot
        if(!is.null(Zr))
            oo <- oo[-(seq_len(zrank))] - ncol(Zr)
        oo <- oo[seq_len(yrank)]
        if(length(oo) < ncol(Yr))
            RDA$alias <- colnames(Yr)[-oo]
        RDA$biplot <- cor(Yr[, oo, drop = FALSE],
                          sol$u[, 1:rank, drop = FALSE])
        RDA$rank <- rank
        RDA$qrank <- Q$rank
        RDA$tot.chi <- sum(RDA$eig)
        RDA$QR <- Q
        RDA$envcentre <- attr(Yr, "scaled:center")
        RDA$Xbar <- Xbar
    } else {
        RDA <- list(eig = 0, rank = rank, qrank = Q$rank, tot.chi = 0,
                    QR = Q, Xbar = Xbar)
        u <- matrix(0, nrow=nrow(sol$u), ncol=0)
        v <- matrix(0, nrow=nrow(sol$v), ncol=0)
        RDA$u <- RDA$wa <- u
        RDA$v <- v
        RDA$biplot <- matrix(0, 0, 0)
        RDA$alias <- colnames(Yr)
    }
    RDA$colsum <- SD
    class(RDA) <- "rdaFit"
    RDA
}

`scores.rdaFit` <- function(x, choices = c(1,2),
                          display = c("sp","wa","cn"), scaling = 2, ...) {
    display <- match.arg(display, c("sites", "species", "wa",
                                    "lc", "bp", "cn"),
                         several.ok = TRUE)
    if("sites" %in% display)
        display[display == "sites"] <- "wa"
    if("species" %in% display)
        display[display == "species"] <- "sp"
    tabula <- c("species", "sites", "constraints", "biplot",
                "centroids")
    names(tabula) <- c("sp", "wa", "lc", "bp", "cn")
    take <- tabula[display]
    slam <- sqrt(x$lambda[choices])
    rnk <- x$rank
    sol <- list()
    if ("species" %in% take) {
        v <- x$v[, choices, drop = FALSE]
        if (scaling) {
            scal <- list(1, slam, sqrt(slam))[[abs(scaling)]]
            v <- sweep(v, 2, scal, "*")
            if (scaling < 0) {
                scal <- sqrt(1/(1 - slam^2))
                v <- sweep(v, 2, scal, "*")
            }
        }
        sol$species <- v
    }
    if ("sites" %in% take) {
        wa <- x$wa[, choices, drop = FALSE]
        if (scaling) {
            scal <- list(slam, 1, sqrt(slam))[[abs(scaling)]]
            wa <- sweep(wa, 2, scal, "*")
            if (scaling < 0) {
                scal <- sqrt(1/(1 - slam^2))
                wa <- sweep(wa, 2, scal, "*")
            }
        }
        sol$sites <- wa
    }
    if ("constraints" %in% take) {
        u <- x$u[, choices, drop = FALSE]
        if (scaling) {
            scal <- list(slam, 1, sqrt(slam))[[abs(scaling)]]
            u <- sweep(u, 2, scal, "*")
            if (scaling < 0) {
                scal <- sqrt(1/(1 - slam^2))
                u <- sweep(u, 2, scal, "*")
            }
        }
        sol$constraints <- u
    }
    if ("biplot" %in% take && !is.null(x$biplot)) {
        b <- matrix(0, nrow(x$biplot), length(choices))
        b[, choices <= rnk] <- x$biplot[, choices[choices <=
            rnk]]
        colnames(b) <- colnames(x$u)[choices]
        rownames(b) <- rownames(x$biplot)
        sol$biplot <- b
    }
    if ("centroids" %in% take) {
        if (is.null(x$centroids))
            sol$centroids <- NA
        else {
            cn <- matrix(0, nrow(x$centroids), length(choices))
            cn[, choices <= rnk] <- x$centroids[, choices[choices <=
                 rnk]]
            colnames(cn) <- colnames(x$u)[choices]
            rownames(cn) <- rownames(x$centroids)
            if (scaling) {
                scal <- list(slam, 1, sqrt(slam))[[abs(scaling)]]
                cn <- sweep(cn, 2, scal, "*")
                if (scaling < 0) {
                    scal <- sqrt(1/(1 - slam^2))
                    cn <- sweep(cn, 2, scal, "*")
                }
            }
            sol$centroids <- cn
        }
    }
    ## Only one type of scores: return a matrix instead of a list
    if (length(sol) == 1)
        sol <- sol[[1]]
    sol
}
