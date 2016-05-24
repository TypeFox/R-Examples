gen.lms <-
function (dat, p = 0.5, bet, epsilon = 0, k.u = 0) 
{
    n <- nrow(dat)
    if (k.u == 0) {
        k.u <- kl.ku(n, p, bet, epsilon)[2]
    }
    else {
        if (k.u < 0) {
            stop("k.u must be a positive integer!\n")
        }
    }
    dat.b <- dat[rowSums(abs(dat)) < Inf, ]
    n.b <- nrow(dat.b)
    dat.b0 <- dat[rowSums(abs(dat[, c(3, 4)])) < Inf, ]
    n.b0 <- nrow(dat.b0)
    if (n.b0 < k.u) {
        stop("Too many unbounded data!\n")
    }
    else {
        dat.unique <- unique(dat, MARGIN = 1)
        m <- nrow(dat.unique)
        ind <- 1:m
        ind.vi <- 0
        for (i in 1:(m - 1)) {
            ind.vi <- c(ind.vi, rep(i, (m - i)))
        }
        ind.vj <- 0
        for (j in 2:m) {
            ind.vj <- c(ind.vj, ind[j:m])
        }
        ind.dat <- matrix(c(ind.vi[-1], ind.vj[-1]), ncol = 2)
        b.pot <- matrix(0, choose(m, 2) + 1, 4)
        for (i in 1:choose(m, 2)) {
            if (dat.unique[ind.dat[i, 1], 1] > dat.unique[ind.dat[i, 
                2], 1]) {
                if (dat.unique[ind.dat[i, 1], 4] > dat.unique[ind.dat[i, 
                  2], 4]) {
                  b.pot[i, 1] <- (dat.unique[ind.dat[i, 1], 4] - 
                    dat.unique[ind.dat[i, 2], 4])/(dat.unique[ind.dat[i, 
                    1], 1] - dat.unique[ind.dat[i, 2], 1])
                }
                if (dat.unique[ind.dat[i, 1], 3] < dat.unique[ind.dat[i, 
                  2], 3]) {
                  b.pot[i, 3] <- (dat.unique[ind.dat[i, 1], 3] - 
                    dat.unique[ind.dat[i, 2], 3])/(dat.unique[ind.dat[i, 
                    1], 1] - dat.unique[ind.dat[i, 2], 1])
                }
            }
            else {
                if (dat.unique[ind.dat[i, 2], 4] > dat.unique[ind.dat[i, 
                  1], 4]) {
                  b.pot[i, 1] <- (dat.unique[ind.dat[i, 2], 4] - 
                    dat.unique[ind.dat[i, 1], 4])/(dat.unique[ind.dat[i, 
                    2], 1] - dat.unique[ind.dat[i, 1], 1])
                }
                if (dat.unique[ind.dat[i, 2], 3] < dat.unique[ind.dat[i, 
                  1], 3]) {
                  b.pot[i, 3] <- (dat.unique[ind.dat[i, 2], 3] - 
                    dat.unique[ind.dat[i, 1], 3])/(dat.unique[ind.dat[i, 
                    2], 1] - dat.unique[ind.dat[i, 1], 1])
                }
            }
            if (dat.unique[ind.dat[i, 1], 2] > dat.unique[ind.dat[i, 
                2], 2]) {
                if (dat.unique[ind.dat[i, 1], 4] < dat.unique[ind.dat[i, 
                  2], 4]) {
                  b.pot[i, 2] <- (dat.unique[ind.dat[i, 1], 4] - 
                    dat.unique[ind.dat[i, 2], 4])/(dat.unique[ind.dat[i, 
                    1], 2] - dat.unique[ind.dat[i, 2], 2])
                }
                if (dat.unique[ind.dat[i, 1], 3] > dat.unique[ind.dat[i, 
                  2], 3]) {
                  b.pot[i, 4] <- (dat.unique[ind.dat[i, 1], 3] - 
                    dat.unique[ind.dat[i, 2], 3])/(dat.unique[ind.dat[i, 
                    1], 2] - dat.unique[ind.dat[i, 2], 2])
                }
            }
            else {
                if (dat.unique[ind.dat[i, 2], 4] < dat.unique[ind.dat[i, 
                  1], 4]) {
                  b.pot[i, 2] <- (dat.unique[ind.dat[i, 2], 4] - 
                    dat.unique[ind.dat[i, 1], 4])/(dat.unique[ind.dat[i, 
                    2], 2] - dat.unique[ind.dat[i, 1], 2])
                }
                if (dat.unique[ind.dat[i, 2], 3] > dat.unique[ind.dat[i, 
                  1], 3]) {
                  b.pot[i, 4] <- (dat.unique[ind.dat[i, 2], 3] - 
                    dat.unique[ind.dat[i, 1], 3])/(dat.unique[ind.dat[i, 
                    1], 2] - dat.unique[ind.dat[i, 1], 2])
                }
            }
        }
        b.pot <- unique(as.vector(round(b.pot, 10)))
        b.pot <- na.omit(b.pot)
        if (max(b.pot) == Inf) {
            b.pot <- b.pot[-which(b.pot == Inf)]
        }
        if (min(b.pot) == -Inf) {
            b.pot <- b.pot[-which(b.pot == -Inf)]
        }
    }
    b.val <- b.pot
    if (n.b < k.u) {
        b.pot <- 0
    }
    a.pot <- rep(0, length(b.pot))
    q.pot <- rep(0, length(b.pot))
    for (i in 1:length(b.pot)) {
        if (b.pot[i] == 0) {
            d <- dat.b0[order(dat.b0[, 4]), c(3, 4)]
            d.u <- d[, 2]
            d.l <- d[, 1]
            d.l <- sort(d.l)
            diff <- matrix(0, (n.b0 - k.u + 1), 3)
            for (k in 1:nrow(diff)) {
                d.u.k <- d[d[, 1] >= d.l[k], 2]
                diff[k, 1] <- d.u.k[k.u] - d.l[k]
                diff[k, 2] <- (d.u.k[k.u] + d.l[k])/2
                diff[k, 3] <- (d.u.k[k.u] - d.l[k])/2
            }
            diff <- unique(diff, MARGIN = 1)
            k.opt <- which(diff[, 1] == min(diff[, 1]))
            if (length(k.opt) > 1) {
                for (l in 2:length(k.opt)) {
                  b.pot <- c(b.pot, b.pot[i])
                  a.pot <- c(a.pot, diff[k.opt[l], 2])
                  q.pot <- c(q.pot, diff[k.opt[l], 3])
                }
            }
            a.pot[i] <- diff[k.opt[1], 2]
            q.pot[i] <- diff[k.opt[1], 3]
        }
        else {
            d <- matrix(0, n.b, 2)
            if (b.pot[i] > 0) {
                d[, 1] <- dat.b[, 3] - b.pot[i] * dat.b[, 2]
                d[, 2] <- dat.b[, 4] - b.pot[i] * dat.b[, 1]
            }
            else {
                d[, 1] <- dat.b[, 3] - b.pot[i] * dat.b[, 1]
                d[, 2] <- dat.b[, 4] - b.pot[i] * dat.b[, 2]
            }
            d <- d[order(d[, 2]), ]
            d.u <- d[, 2]
            d.l <- d[, 1]
            d.l <- sort(d.l)
            diff <- matrix(0, (n.b - k.u + 1), 3)
            for (k in 1:nrow(diff)) {
                d.u.k <- d[d[, 1] >= d.l[k], 2]
                diff[k, 1] <- d.u.k[k.u] - d.l[k]
                diff[k, 2] <- (d.u.k[k.u] + d.l[k])/2
                diff[k, 3] <- (d.u.k[k.u] - d.l[k])/2
            }
            diff <- unique(diff, MARGIN = 1)
            k.opt <- which(diff[, 1] == min(diff[, 1]))
            if (length(k.opt) > 1) {
                for (l in 2:length(k.opt)) {
                  b.pot <- c(b.pot, b.pot[i])
                  a.pot <- c(a.pot, diff[k.opt[l], 2])
                  q.pot <- c(q.pot, diff[k.opt[l], 3])
                }
            }
            a.pot[i] <- diff[k.opt[1], 2]
            q.pot[i] <- diff[k.opt[1], 3]
        }
    }
    preresult <- cbind(a.pot, b.pot, q.pot)
    preresult <- unique(preresult)
    attr(preresult, "dimnames")[[2]] <- c("a.lrm", "b.lrm", "q.lrm")
    list(lrm = preresult[preresult[, 3] == min(preresult[, 3]), 
        ], b.val = sort(b.val), dat.unique = dat.unique, ind.dat = ind.dat)
}
idf.create <-
function (dat, var.labels = NULL) 
{
    if (!(is.data.frame(dat))) 
        stop("The data must be provided as a data frame! \n")
    m <- ncol(dat)/2
    if (ceiling(m) != m) 
        stop("The data frame must contain two neighboring columns (with lower and upper endpoints, respectively) for each variable! \n")
    idf <- vector("list", m)
    var.lab <- NULL
    for (i in 1:m) {
        idf[[i]] <- dat[, c((2 * i - 1), 2 * i)]
    }
    if (!is.null(var.labels)) {
        for (i in 1:m) {
            attr(idf[[i]], "names") <- c(paste(var.labels[i], 
                "l", sep = "."), paste(var.labels[i], "u", sep = "."))
        }
        attr(idf, "names") <- var.labels
    }
    else {
        for (i in 1:m) {
            attr(idf[[i]], "names") <- c(paste(paste("var", i, 
                sep = ""), "l", sep = "."), paste(paste("var", 
                i, sep = ""), "u", sep = "."))
            var.lab <- c(var.lab, paste("var", i, sep = ""))
        }
        attr(idf, "names") <- var.lab
    }
    idf$n <- nrow(dat)
    class(idf) <- "idf"
    idf
}
kl.ku <-
function (n, p = 0.5, bet, epsilon = 0) 
{
    if ((max(p, (1 - p)) + epsilon)^n > bet) {
        stop("k.l and k.u are not defined ! \n")
    }
    lambda <- function(s, t) {
        (s/t)^(-s) * ((1 - s)/(1 - t))^(s - 1)
    }
    k.l <- Find(function(k) {
        lambda(k/n, p - epsilon) <= bet^(1/n)
    }, 0:floor((p - epsilon) * n), right = T)
    k.u <- Find(function(k) {
        lambda(k/n, p + epsilon) <= bet^(1/n)
    }, (floor((p + epsilon) * n) + 1):n)
    c(k.l, k.u)
}
plot.idf <-
function (x, y = NULL, ..., var = NULL, typ = "hist", k.x = 1, 
    k.y = 1, inf.margin = 10, p.cex = 1, col.lev = 15, plot.grid = FALSE, 
    x.adj = 0.5, x.padj = 3, y.las = 0, y.adj = 1, y.padj = 0, 
    x.lim = c(0, 0), y.lim = c(0, 0), x.lab = "X", y.lab = "Y") 
{
    dat.idf <- x
    if (!is.null(var)) {
        dat <- cbind(dat.idf[[which(names(dat.idf) == var[1])]], 
            dat.idf[[which(names(dat.idf) == var[2])]])
    }
    else {
        dat <- cbind(dat.idf[[1]], dat.idf[[2]])
    }
    if (x.lim[1] == 0 & x.lim[2] == 0) {
        x.min <- floor(min(c(dat[dat[, 1] != -Inf, 1], dat[dat[, 
            2] != Inf, 2])))
        x.max <- ceiling(max(c(dat[dat[, 1] != -Inf, 1], dat[dat[, 
            2] != Inf, 2])))
    }
    else {
        if (min(c(dat[dat[, 1] != -Inf, 1], dat[dat[, 2] != Inf, 
            2])) < x.lim[1] | max(c(dat[dat[, 1] != -Inf, 1], 
            dat[dat[, 2] != Inf, 2])) > x.lim[2]) {
            stop("There are observations outside the chosen x-limits! \n")
        }
        x.min <- x.lim[1]
        x.max <- x.lim[2]
    }
    if (y.lim[1] == 0 & y.lim[2] == 0) {
        y.min <- floor(min(c(dat[dat[, 3] != -Inf, 3], dat[dat[, 
            4] != Inf, 4])))
        y.max <- ceiling(max(c(dat[dat[, 3] != -Inf, 3], dat[dat[, 
            4] != Inf, 4])))
    }
    else {
        if (min(c(dat[dat[, 3] != -Inf, 3], dat[dat[, 4] != Inf, 
            4])) < y.lim[1] | max(c(dat[dat[, 3] != -Inf, 3], 
            dat[dat[, 4] != Inf, 4])) > y.lim[2]) {
            stop("There are observations outside the chosen y-limits! \n")
        }
        y.min <- y.lim[1]
        y.max <- y.lim[2]
    }
    n <- dat.idf$n
    x.range <- x.max - x.min
    y.range <- y.max - y.min
    X.min <- round((x.min - 0.04 * x.range) * k.x)/k.x
    X.max <- round((x.max + 0.04 * x.range) * k.x)/k.x
    Y.min <- round((y.min - 0.04 * y.range) * k.y)/k.y
    Y.max <- round((y.max + 0.04 * y.range) * k.y)/k.y
    if (typ == "draft") {
        dat[dat[, 1] == -Inf, 1] <- X.min - inf.margin/k.x
        dat[dat[, 2] == Inf, 2] <- X.max + inf.margin/k.x
        dat[dat[, 3] == -Inf, 3] <- Y.min - inf.margin/k.y
        dat[dat[, 4] == Inf, 4] <- Y.max + inf.margin/k.y
        dat[, 1] <- round(dat[, 1] * k.x)/k.x
        dat[, 2] <- round(dat[, 2] * k.x)/k.x
        dat[, 3] <- round(dat[, 3] * k.y)/k.y
        dat[, 4] <- round(dat[, 4] * k.y)/k.y
        dat <- dat[order(dat[, 1]), ]
        x.steps <- rep(0, n)
        y.steps <- rep(0, n)
        for (i in 1:nrow(dat)) {
            x.steps[i] <- round((dat[i, 2] - dat[i, 1]) * k.x, 
                10)
            y.steps[i] <- round((dat[i, 4] - dat[i, 3]) * k.y, 
                10)
        }
        plot(mean(dat[, 1]), mean(dat[, 3]), type = "p", pch = 15, 
            col = 0, xlim = c(x.min, x.max), ylim = c(y.min, 
                y.max), xlab = " ", ylab = " ", las = y.las)
        mtext(x.lab, side = 1, las = 1, adj = x.adj, padj = x.padj)
        mtext(y.lab, side = 2, las = y.las, adj = y.adj, padj = y.padj)
        for (i in 1:nrow(dat)) {
            if (x.steps[i] > 0 & y.steps[i] > 0) {
                Z.i <- matrix(0, nrow = (x.steps[i] * y.steps[i]), 
                  2)
                if (x.steps[i] == 1) {
                  Z.i[, 1] <- rep(dat[i, 1] + 1/(2 * k.x), times = y.steps[i])
                }
                else {
                  Z.i[, 1] <- rep(seq(dat[i, 1] + 1/(2 * k.x), 
                    dat[i, 2] - 1/(2 * k.x), 1/k.x), times = y.steps[i])
                }
                if (y.steps[i] == 1) {
                  Z.i[, 2] <- rep(dat[i, 3] + 1/(2 * k.y), each = x.steps[i])
                }
                else {
                  Z.i[, 2] <- rep(seq(dat[i, 3] + 1/(2 * k.y), 
                    dat[i, 4] - 1/(2 * k.y), 1/k.y), each = x.steps[i])
                }
                points(Z.i[, 1], Z.i[, 2], pch = 15, col = "lightgray", 
                  cex = p.cex)
            }
            points(dat[i, 2], dat[i, 4], pch = 20, cex = 0.25)
            points(dat[i, 2], dat[i, 3], pch = 20, cex = 0.25)
            points(dat[i, 1], dat[i, 4], pch = 20, cex = 0.25)
            points(dat[i, 1], dat[i, 3], pch = 20, cex = 0.25)
            segments(dat[i, 1], dat[i, 4], dat[i, 2], dat[i, 
                4], lty = 1, lwd = 2)
            segments(dat[i, 1], dat[i, 3], dat[i, 2], dat[i, 
                3], lty = 1, lwd = 2)
            segments(dat[i, 1], dat[i, 3], dat[i, 1], dat[i, 
                4], lty = 1, lwd = 2)
            segments(dat[i, 2], dat[i, 3], dat[i, 2], dat[i, 
                4], lty = 1, lwd = 2)
        }
    }
    else {
        x.plot.limits <- sort(unique(round(c(dat[, 1], dat[, 
            2]) * k.x)/k.x))
        y.plot.limits <- sort(unique(round(c(dat[, 3], dat[, 
            4]) * k.y)/k.y))
        dat[dat[, 1] == -Inf, 1] <- X.min
        dat[dat[, 2] == Inf, 2] <- X.max
        dat[dat[, 3] == -Inf, 3] <- Y.min
        dat[dat[, 4] == Inf, 4] <- Y.max
        dat[, 1] <- round(dat[, 1] * k.x)/k.x
        dat[, 2] <- round(dat[, 2] * k.x)/k.x
        dat[, 3] <- round(dat[, 3] * k.y)/k.y
        dat[, 4] <- round(dat[, 4] * k.y)/k.y
        dat <- dat[order(dat[, 1]), ]
        x.steps <- rep(0, n)
        y.steps <- rep(0, n)
        for (i in 1:nrow(dat)) {
            x.steps[i] <- round((dat[i, 2] - dat[i, 1]) * k.x, 
                10)
            y.steps[i] <- round((dat[i, 4] - dat[i, 3]) * k.y, 
                10)
        }
        x.steps[which(x.steps == 0)] <- 1
        y.steps[which(y.steps == 0)] <- 1
        Z.i <- matrix(0, x.steps[1] * y.steps[1], 2)
        Z.i[, 1] <- rep(seq(round(dat[1, 1] + 1/(2 * k.x), 10), 
            round(dat[1, 1] + (x.steps[1]) * 1/k.x - 1/(2 * k.x), 
                10), round(1/k.x, 10)), times = y.steps[1])
        Z.i[, 2] <- rep(seq(round(dat[1, 3] + 1/(2 * k.y), 10), 
            round(dat[1, 3] + (y.steps[1]) * 1/k.y - 1/(2 * k.y), 
                10), round(1/k.y, 10)), each = x.steps[1])
        Z <- Z.i
        for (i in 2:nrow(dat)) {
            Z.i <- matrix(0, nrow = (x.steps[i] * y.steps[i]), 
                2)
            Z.i[, 1] <- rep(seq(round(dat[i, 1] + 1/(2 * k.x), 
                10), round(dat[i, 1] + (x.steps[i]) * 1/k.x - 
                1/(2 * k.x), 10), round(1/k.x, 10)), times = y.steps[i])
            Z.i[, 2] <- rep(seq(round(dat[i, 3] + 1/(2 * k.y), 
                10), round(dat[i, 3] + (y.steps[i]) * 1/k.y - 
                1/(2 * k.y), 10), round(1/k.y, 10)), each = x.steps[i])
            Z <- rbind(Z, Z.i)
        }
        X <- round(seq(X.min - 1/(2 * k.x), X.max + 1/(2 * k.x), 
            1/k.x), 10)
        Y <- round(seq(Y.min - 1/(2 * k.y), Y.max + 1/(2 * k.y), 
            1/k.y), 10)
        Z.image <- cbind(rep(X, times = length(Y)), rep(Y, each = length(X)))
        Z <- rbind(Z, Z.image)
        tab <- ftable(round(Z[, 1], 10), round(Z[, 2], 10))
        max.tab <- max(tab)
        image(X, Y, tab, col = gray((col.lev:0)/col.lev), xlim = c(X.min, 
            X.max), ylim = c(Y.min, Y.max), zlim = c(1, max.tab), 
            xlab = " ", ylab = " ", las = y.las)
        mtext(x.lab, side = 1, las = 1, adj = x.adj, padj = x.padj)
        mtext(y.lab, side = 2, las = y.las, adj = y.adj, padj = y.padj)
        if (plot.grid == TRUE) {
            abline(v = x.plot.limits, lty = 2)
            abline(h = y.plot.limits, lty = 2)
        }
        abline(v = c(X.min, X.max))
        abline(h = c(Y.min, Y.max))
    }
}
plot.s.linlir <-
function (x, y = NULL, ..., typ, para.typ = "polygon", b.grid = 500, 
    nb.func = 1000, seed.func = NULL, pl.lrm = TRUE, pl.band = FALSE, 
    lrm.col = "blue", pl.dat = FALSE, pl.dat.typ = "hist", k.x = 1, 
    k.y = 1, inf.margin = 10, p.cex = 1, col.lev = 15, plot.grid = FALSE, 
    x.adj = 0.5, x.padj = 3, y.las = 0, y.adj = 1, y.padj = 0, 
    x.lim = c(0, 0), y.lim = c(0, 0), x.lab = " ", y.lab = " ") 
{
    x.s.linlir <- x
    if (typ == "para") {
        if (sum(abs(x.s.linlir$a.undom)) == Inf | sum(abs(x.s.linlir$b.undom)) == 
            Inf) {
            if ((x.lim[1] == 0 & x.lim[2] == 0) | (y.lim[1] == 
                0 & y.lim[2] == 0)) {
                stop("Parameter set unbounded. Choose plot limits x.lim and y.lim !\n")
            }
        }
        if (x.lim[1] == 0 & x.lim[2] == 0) {
            x.min <- floor(x.s.linlir$b.undom[[1]])
            x.max <- ceiling(x.s.linlir$b.undom[[2]])
        }
        else {
            x.min <- x.lim[1]
            x.max <- x.lim[2]
        }
        if (y.lim[1] == 0 & y.lim[2] == 0) {
            y.min <- floor(x.s.linlir$a.undom[[1]])
            y.max <- ceiling(x.s.linlir$a.undom[[2]])
        }
        else {
            y.min <- y.lim[1]
            y.max <- y.lim[2]
        }
        if (x.lab == " ") {
            x.lab <- "b"
        }
        if (y.lab == " ") {
            y.lab <- "a"
        }
        if (para.typ == "polygon") {
            b.d <- inf.margin/100 * (x.max - x.min)
            b.pot <- seq(x.min - b.d, x.max + b.d, by = (x.max - 
                x.min + 2 * b.d)/(b.grid - 1))
            n <- x.s.linlir$n
            k.l <- x.s.linlir$config$k.l
            a.l.plot <- matrix(NA, nrow = (n - k.l), ncol = length(b.pot))
            a.u.plot <- matrix(NA, nrow = (n - k.l), ncol = length(b.pot))
            for (i in 1:length(b.pot)) {
                a.int <- undom.a(x.s.linlir$dat, b = b.pot[i], 
                  x.s.linlir$q.lrm, x.s.linlir$config$p, x.s.linlir$config$bet, 
                  x.s.linlir$config$epsilon)
                a.l.plot[, i] <- a.int[[1]][, 1]
                a.u.plot[, i] <- a.int[[1]][, 2]
            }
            Y.min <- y.min - inf.margin/100 * (y.max - y.min)
            Y.max <- y.max + inf.margin/100 * (y.max - y.min)
            a.l.plot[which(a.l.plot <= Y.min)] <- Y.min
            a.u.plot[which(a.u.plot >= Y.max)] <- Y.max
            plot(x.min, y.max, type = "n", xlim = c(x.min, x.max), 
                ylim = c(y.min, y.max), xlab = " ", ylab = " ", 
                las = y.las)
            mtext(x.lab, side = 1, las = 1, adj = x.adj, padj = x.padj)
            mtext(y.lab, side = 2, las = y.las, adj = y.adj, 
                padj = y.padj)
            for (j in 1:(n - k.l)) {
                polygon(c(b.pot[a.l.plot[j, ] <= a.u.plot[j, 
                  ]], rev(b.pot[a.l.plot[j, ] <= a.u.plot[j, 
                  ]])), c(a.l.plot[j, a.l.plot[j, ] <= a.u.plot[j, 
                  ]], rev(a.u.plot[j, a.l.plot[j, ] <= a.u.plot[j, 
                  ]])), col = "darkgrey", lty = 0)
            }
            if (!is.matrix(x.s.linlir$undom.para)) {
                points(x.s.linlir$undom.para[2], x.s.linlir$undom.para[1], 
                  pch = 19, cex = 0.5, col = "darkgrey")
            }
        }
        else {
            undom.para <- x.s.linlir$undom.para
            if (is.matrix(undom.para)) {
                plot(undom.para[, 2], undom.para[, 1], pch = 19, 
                  cex = 0.5, col = "darkgrey", xlim = c(x.min, 
                    x.max), ylim = c(y.min, y.max), xlab = " ", 
                  ylab = " ", las = y.las)
                mtext(x.lab, side = 1, las = 1, adj = x.adj, 
                  padj = x.padj)
                mtext(y.lab, side = 2, las = y.las, adj = y.adj, 
                  padj = y.padj)
            }
            else {
                plot(undom.para[2], undom.para[1], pch = 19, 
                  cex = 0.5, col = "darkgrey", xlim = c(x.min, 
                    x.max), ylim = c(y.min, y.max), xlab = " ", 
                  ylab = " ", las = y.las)
                mtext(x.lab, side = 1, las = 1, adj = x.adj, 
                  padj = x.padj)
                mtext(y.lab, side = 2, las = y.las, adj = y.adj, 
                  padj = y.padj)
            }
        }
        if (pl.lrm == TRUE) {
            if (!is.vector(x.s.linlir$f.lrm)) {
                for (j in 1:nrow(x.s.linlir$f.lrm)) {
                  points(x.s.linlir$f.lrm[j, 2], x.s.linlir$f.lrm[j, 
                    1], pch = 19, col = lrm.col, cex = 1)
                }
                print("f.lrm is not unique !\n")
            }
            else {
                points(x.s.linlir$f.lrm[2], x.s.linlir$f.lrm[1], 
                  pch = 19, col = lrm.col, cex = 1)
            }
        }
    }
    else {
        dat <- x.s.linlir$dat
        if (x.lim[1] == 0 & x.lim[2] == 0) {
            x.min <- floor(min(c(dat[dat[, 1] != -Inf, 1], dat[dat[, 
                2] != Inf, 2])))
            x.max <- ceiling(max(c(dat[dat[, 1] != -Inf, 1], 
                dat[dat[, 2] != Inf, 2])))
        }
        else {
            x.min <- x.lim[1]
            x.max <- x.lim[2]
        }
        if (y.lim[1] == 0 & y.lim[2] == 0) {
            y.min <- floor(min(c(dat[dat[, 3] != -Inf, 3], dat[dat[, 
                4] != Inf, 4])))
            y.max <- ceiling(max(c(dat[dat[, 3] != -Inf, 3], 
                dat[dat[, 4] != Inf, 4])))
        }
        else {
            y.min <- y.lim[1]
            y.max <- y.lim[2]
        }
        x.d <- inf.margin/100 * (x.max - x.min)
        if (x.lab == " ") {
            x.lab <- "X"
        }
        if (y.lab == " ") {
            y.lab <- "Y"
        }
        if (typ == "lrm") {
            if (!is.vector(x.s.linlir$f.lrm)) {
                if (pl.dat == TRUE) {
                  dat.idf <- idf.create(dat)
                  plot(dat.idf, typ = pl.dat.typ, k.x = k.x, 
                    k.y = k.y, inf.margin = inf.margin, p.cex = p.cex, 
                    col.lev = col.lev, plot.grid = plot.grid, 
                    x.adj = x.adj, x.padj = x.padj, y.las = y.las, 
                    y.adj = y.adj, y.padj = y.padj, x.lim = c(x.min, 
                      x.max), y.lim = c(y.min, y.max), x.lab = x.lab, 
                    y.lab = y.lab)
                  for (j in 1:nrow(x.s.linlir$f.lrm)) {
                    curve(x.s.linlir$f.lrm[j, 1] + x.s.linlir$f.lrm[j, 
                      2] * x, x.min - x.d, x.max + x.d, add = T, 
                      lty = 1, col = lrm.col, lwd = 2)
                  }
                }
                else {
                  plot(min(dat[, 1]), max(dat[, 4]), type = "n", 
                    xlim = c(x.min, x.max), ylim = c(y.min, y.max), 
                    xlab = " ", ylab = " ", las = y.las)
                  mtext(x.lab, side = 1, las = 1, adj = x.adj, 
                    padj = x.padj)
                  mtext(y.lab, side = 2, las = y.las, adj = y.adj, 
                    padj = y.padj)
                  for (j in 1:nrow(x.s.linlir$f.lrm)) {
                    curve(x.s.linlir$f.lrm[j, 1] + x.s.linlir$f.lrm[j, 
                      2] * x, x.min - x.d, x.max + x.d, add = T, 
                      lty = 1, col = lrm.col, lwd = 2)
                  }
                }
                print("f.lrm is not unique !\n")
            }
            else {
                if (pl.dat == TRUE) {
                  dat.idf <- idf.create(dat)
                  plot(dat.idf, typ = pl.dat.typ, k.x = k.x, 
                    k.y = k.y, inf.margin = inf.margin, p.cex = p.cex, 
                    col.lev = col.lev, plot.grid = plot.grid, 
                    x.adj = x.adj, x.padj = x.padj, y.las = y.las, 
                    y.adj = y.adj, y.padj = y.padj, x.lim = c(x.min, 
                      x.max), y.lim = c(y.min, y.max), x.lab = x.lab, 
                    y.lab = y.lab)
                  curve(x.s.linlir$f.lrm[1] + x.s.linlir$f.lrm[2] * 
                    x, x.min - x.d, x.max + x.d, add = T, lty = 1, 
                    col = lrm.col, lwd = 2)
                }
                else {
                  plot(min(dat[, 1]), max(dat[, 4]), type = "n", 
                    xlim = c(x.min, x.max), ylim = c(y.min, y.max), 
                    xlab = " ", ylab = " ", las = y.las)
                  mtext(x.lab, side = 1, las = 1, adj = x.adj, 
                    padj = x.padj)
                  mtext(y.lab, side = 2, las = y.las, adj = y.adj, 
                    padj = y.padj)
                  curve(x.s.linlir$f.lrm[1] + x.s.linlir$f.lrm[2] * 
                    x, x.min - x.d, x.max + x.d, add = T, lty = 1, 
                    col = lrm.col, lwd = 2)
                }
            }
        }
        else {
            if (nb.func > nrow(x.s.linlir$undom.para)) {
                nb.func <- nrow(x.s.linlir$undom.para)
            }
            if (!is.null(seed.func)) {
                set.seed(seed.func)
            }
            undom <- x.s.linlir$undom.para[sample(1:nrow(x.s.linlir$undom.para), 
                nb.func), ]
            if (pl.dat == TRUE) {
                dat.idf <- idf.create(dat)
                plot(dat.idf, typ = pl.dat.typ, k.x = k.x, k.y = k.y, 
                  inf.margin = inf.margin, p.cex = p.cex, col.lev = col.lev, 
                  plot.grid = plot.grid, x.adj = x.adj, x.padj = x.padj, 
                  y.las = y.las, y.adj = y.adj, y.padj = y.padj, 
                  x.lim = c(x.min, x.max), y.lim = c(y.min, y.max), 
                  x.lab = x.lab, y.lab = y.lab)
                for (i in 1:nrow(undom)) {
                  curve(undom[i, 1] + undom[i, 2] * x, x.min - 
                    x.d, x.max + x.d, add = T, lty = 1, col = "darkgrey")
                }
            }
            else {
                plot(min(dat[, 1]), max(dat[, 4]), type = "n", 
                  xlim = c(x.min, x.max), ylim = c(y.min, y.max), 
                  xlab = " ", ylab = " ", las = y.las)
                mtext(x.lab, side = 1, las = 1, adj = x.adj, 
                  padj = x.padj)
                mtext(y.lab, side = 2, las = y.las, adj = y.adj, 
                  padj = y.padj)
                for (i in 1:nrow(undom)) {
                  curve(undom[i, 1] + undom[i, 2] * x, x.min - 
                    x.d, x.max + x.d, add = T, lty = 1, col = "darkgrey")
                }
            }
            if (pl.lrm == TRUE) {
                if (!is.vector(x.s.linlir$f.lrm)) {
                  for (j in 1:nrow(x.s.linlir$f.lrm)) {
                    curve(x.s.linlir$f.lrm[j, 1] + x.s.linlir$f.lrm[j, 
                      2] * x, x.min - x.d, x.max + x.d, add = T, 
                      lty = 1, col = lrm.col, lwd = 2)
                  }
                  print("f.lrm is not unique !\n")
                }
                else {
                  curve(x.s.linlir$f.lrm[1] + x.s.linlir$f.lrm[2] * 
                    x, x.min - x.d, x.max + x.d, add = T, lty = 1, 
                    col = lrm.col, lwd = 2)
                }
            }
        }
        if (pl.band == TRUE) {
            if (!is.vector(x.s.linlir$f.lrm)) {
                for (j in 1:nrow(x.s.linlir$f.lrm)) {
                  curve(x.s.linlir$f.lrm[j, 1] + x.s.linlir$q.lrm + 
                    x.s.linlir$f.lrm[j, 2] * x, x.min - x.d, 
                    x.max + x.d, add = T, lty = 2, col = lrm.col, 
                    lwd = 2)
                  curve(x.s.linlir$f.lrm[j, 1] + x.s.linlir$q.lrm + 
                    x.s.linlir$f.lrm[j, 2] * x, x.min - x.d, 
                    x.max + x.d, add = T, lty = 2, col = lrm.col, 
                    lwd = 2)
                }
                print("f.lrm is not unique !\n")
            }
            else {
                curve(x.s.linlir$f.lrm[1] + x.s.linlir$q.lrm + 
                  x.s.linlir$f.lrm[2] * x, x.min - x.d, x.max + 
                  x.d, add = T, lty = 2, col = lrm.col, lwd = 2)
                curve(x.s.linlir$f.lrm[1] - x.s.linlir$q.lrm + 
                  x.s.linlir$f.lrm[2] * x, x.min - x.d, x.max + 
                  x.d, add = T, lty = 2, col = lrm.col, lwd = 2)
            }
        }
    }
}
print.s.linlir <-
function (x, ...) 
{
    x.s.linlir <- x
    cat("\nSimple linear LIR analysis\n")
    cat("\nCall:\n")
    print(x.s.linlir$call)
    cat("\nLIR settings:\n")
    cat(paste("p:", eval(x.s.linlir$config$p), "   beta:", eval(x.s.linlir$config$bet), 
        "   epsilon:", eval(x.s.linlir$config$epsilon), "   k.l:", 
        eval(x.s.linlir$config$k.l), "   k.u:", eval(x.s.linlir$config$k.u), 
        " \n"))
    cat(paste("confidence level of each confidence interval:", 
        eval(round(x.s.linlir$config$conf.lev.ci * 100, 2)), 
        "% \n"))
}
s.linlir <-
function (dat.idf, var = NULL, p = 0.5, bet, epsilon = 0, a.grid = 100) 
{
    if (class(dat.idf) != "idf") {
        stop("The data must be provided as an *idf* object !\n")
    }
    if (!is.null(var)) {
        dat <- cbind(dat.idf[[which(names(dat.idf) == var[1])]], 
            dat.idf[[which(names(dat.idf) == var[2])]])
    }
    else {
        dat <- cbind(dat.idf[[1]], dat.idf[[2]])
    }
    lrm <- gen.lms(dat, p, bet, epsilon, k.u = 0)
    x.s.linlir <- vector("list", 2)
    if (!is.vector(lrm[[1]])) {
        x.s.linlir[[1]] <- lrm[[1]][, 1:2]
        attr(x.s.linlir, "names")[1] <- "f.lrm"
        x.s.linlir[[2]] <- lrm[[1]][1, 3]
        attr(x.s.linlir, "names")[2] <- "q.lrm"
    }
    else {
        x.s.linlir[[1]] <- lrm[[1]][1:2]
        attr(x.s.linlir, "names")[1] <- "f.lrm"
        x.s.linlir[[2]] <- lrm[[1]][3]
        attr(x.s.linlir, "names")[2] <- "q.lrm"
    }
    b.search <- lrm[[2]]
    dat.unique <- lrm[[3]]
    m <- nrow(dat.unique)
    ind.dat <- lrm[[4]]
    b.inter <- matrix(0, choose(m, 2) + 1, 4)
    for (i in 1:choose(m, 2)) {
        if (dat.unique[ind.dat[i, 1], 1] > dat.unique[ind.dat[i, 
            2], 2]) {
            if (dat.unique[ind.dat[i, 1], 3] < (dat.unique[ind.dat[i, 
                2], 4] + 2 * x.s.linlir$q.lrm)) {
                b.inter[i, 1] <- (dat.unique[ind.dat[i, 1], 3] - 
                  (dat.unique[ind.dat[i, 2], 4] + 2 * x.s.linlir$q.lrm))/(dat.unique[ind.dat[i, 
                  1], 1] - dat.unique[ind.dat[i, 2], 2])
            }
            if ((dat.unique[ind.dat[i, 1], 4] + 2 * x.s.linlir$q.lrm) > 
                dat.unique[ind.dat[i, 2], 3]) {
                b.inter[i, 2] <- ((dat.unique[ind.dat[i, 1], 
                  4] + 2 * x.s.linlir$q.lrm) - dat.unique[ind.dat[i, 
                  2], 3])/(dat.unique[ind.dat[i, 1], 1] - dat.unique[ind.dat[i, 
                  2], 2])
            }
        }
        else {
            if ((dat.unique[ind.dat[i, 2], 4] + 2 * x.s.linlir$q.lrm) < 
                dat.unique[ind.dat[i, 1], 3]) {
                b.inter[i, 1] <- ((dat.unique[ind.dat[i, 2], 
                  4] + 2 * x.s.linlir$q.lrm) - dat.unique[ind.dat[i, 
                  1], 3])/(dat.unique[ind.dat[i, 2], 2] - dat.unique[ind.dat[i, 
                  1], 1])
            }
            if (dat.unique[ind.dat[i, 2], 3] > (dat.unique[ind.dat[i, 
                1], 4] + 2 * x.s.linlir$q.lrm)) {
                b.inter[i, 2] <- (dat.unique[ind.dat[i, 2], 3] - 
                  (dat.unique[ind.dat[i, 1], 4] + 2 * x.s.linlir$q.lrm))/(dat.unique[ind.dat[i, 
                  2], 2] - dat.unique[ind.dat[i, 1], 1])
            }
        }
        if (dat.unique[ind.dat[i, 1], 2] > dat.unique[ind.dat[i, 
            2], 1]) {
            if ((dat.unique[ind.dat[i, 1], 4] + 2 * x.s.linlir$q.lrm) < 
                dat.unique[ind.dat[i, 2], 3]) {
                b.inter[i, 3] <- ((dat.unique[ind.dat[i, 1], 
                  4] + 2 * x.s.linlir$q.lrm) - dat.unique[ind.dat[i, 
                  2], 3])/(dat.unique[ind.dat[i, 1], 2] - dat.unique[ind.dat[i, 
                  2], 1])
            }
            if (dat.unique[ind.dat[i, 1], 3] > (dat.unique[ind.dat[i, 
                2], 4] + 2 * x.s.linlir$q.lrm)) {
                b.inter[i, 4] <- (dat.unique[ind.dat[i, 1], 3] - 
                  (dat.unique[ind.dat[i, 2], 4] + 2 * x.s.linlir$q.lrm))/(dat.unique[ind.dat[i, 
                  1], 2] - dat.unique[ind.dat[i, 2], 1])
            }
        }
        else {
            if (dat.unique[ind.dat[i, 2], 3] < (dat.unique[ind.dat[i, 
                1], 4] + 2 * x.s.linlir$q.lrm)) {
                b.inter[i, 3] <- (dat.unique[ind.dat[i, 2], 3] - 
                  (dat.unique[ind.dat[i, 1], 4] + 2 * x.s.linlir$q.lrm))/(dat.unique[ind.dat[i, 
                  2], 1] - dat.unique[ind.dat[i, 1], 2])
            }
            if ((dat.unique[ind.dat[i, 2], 4] + 2 * x.s.linlir$q.lrm) > 
                dat.unique[ind.dat[i, 1], 3]) {
                b.inter[i, 4] <- ((dat.unique[ind.dat[i, 2], 
                  4] + 2 * x.s.linlir$q.lrm) - dat.unique[ind.dat[i, 
                  1], 3])/(dat.unique[ind.dat[i, 2], 1] - dat.unique[ind.dat[i, 
                  1], 2])
            }
        }
    }
    b.inter <- unique(as.vector(round(b.inter, 10)))
    b.inter <- na.omit(b.inter)
    if (max(b.inter) == Inf) {
        b.inter <- b.inter[-which(b.inter == Inf)]
    }
    if (min(b.inter) == -Inf) {
        b.inter <- b.inter[-which(b.inter == -Inf)]
    }
    b.search <- sort(unique(c(b.search, b.inter)))
    b.length <- length(b.search)
    i <- 1
    while (i < b.length) {
        b.range <- c(b.search[i] - c(100, 1, 0.01), b.search[i], 
            (b.search[i] + b.search[i + 1])/2)
        undom.l <- rep(0, length(b.range))
        for (j in 1:length(b.range)) {
            par.lj <- undom.a(dat, b.range[j], x.s.linlir$q.lrm, 
                p, bet, epsilon)
            undom.l[j] <- as.numeric(is.matrix(par.lj[[2]]))
        }
        if (sum(undom.l) > 0) {
            ind.bl <- i
            i <- b.length
        }
        else {
            i <- i + 1
        }
    }
    if (ind.bl == 1 & undom.l[1] == 1) {
        b.l <- -1e+09
    }
    else {
        b.l <- ceiling(b.search[ind.bl] * 1e+09)/1e+09
    }
    i <- 1
    while (i < b.length) {
        b.range <- c((b.search[b.length - i + 1] + b.search[b.length - 
            i])/2, b.search[b.length - i + 1], b.search[b.length - 
            i + 1] + c(0.01, 1, 100))
        undom.u <- rep(0, length(b.range))
        for (j in 1:length(b.range)) {
            par.uj <- undom.a(dat, b.range[j], x.s.linlir$q.lrm, 
                p, bet, epsilon)
            undom.u[j] <- as.numeric(is.matrix(par.uj[[2]]))
        }
        if (sum(undom.u) > 0) {
            ind.bu <- b.length - i + 1
            i <- b.length
        }
        else {
            i <- i + 1
        }
    }
    if (ind.bu == b.length & undom.u[length(b.range)] == 1) {
        b.u <- 1e+09
    }
    else {
        b.u <- floor(b.search[ind.bu] * 1e+09)/1e+09
    }
    b.step <- (b.u - b.l)/(a.grid * 10 - 1)
    b.range <- c(seq(b.l - b.step, b.u + b.step, b.step), b.l + 
        1e-16, b.u - 1e-16, ceiling(b.l * 1e+06)/1e+06, ceiling(b.l * 
        1000)/1000, floor(b.u * 1e+06)/1e+06, floor(b.u * 1000)/1000)
    para <- undom.para(dat, b.range, a.grid, x.s.linlir$q.lrm, 
        p, bet, epsilon)
    x.s.linlir$a.undom <- round(para[[1]], 10)
    x.s.linlir$b.undom <- round(para[[2]], 10)
    x.s.linlir$undom.para <- para[[3]]
    klku <- kl.ku(dat.idf$n, p, bet, epsilon)
    if (epsilon <= 0) {
        conf.l <- sum(dbinom((klku[1] + 1):klku[2], dat.idf$n, 
            p))
    }
    else {
        if (p <= 0.5) {
            conf.l <- sum(dbinom((klku[1] + 1):klku[2], dat.idf$n, 
                (p + epsilon)))
        }
        else {
            conf.l <- sum(dbinom((klku[1] + 1):klku[2], dat.idf$n, 
                (p - epsilon)))
        }
    }
    as.conf <- pchisq(q = (-2 * log(bet)), df = 1)
    x.s.linlir$config <- list(p = p, beta = bet, epsilon = epsilon, 
        k.l = klku[1], k.u = klku[2], conf.lev.ci = conf.l, as.conf.lev.ci = as.conf)
    x.s.linlir$dat <- dat
    x.s.linlir$n <- dat.idf$n
    x.s.linlir$call <- match.call()
    class(x.s.linlir) <- "s.linlir"
    x.s.linlir
}
summary.idf <-
function (object, ...) 
{
    dat.idf <- object
    cat("\nSummary of interval data frame\n\n")
    cat(dat.idf$n, "observations\n")
    cat((length(dat.idf) - 1), "interval-valued variables \n")
    for (i in 1:(length(dat.idf) - 1)) {
        cat(paste("\nVariable", names(dat.idf)[i], ": \n\n", 
            sep = " "))
        print(summary(dat.idf[[i]]))
    }
}
summary.s.linlir <-
function (object, ...) 
{
    x.s.linlir <- object
    cat("\nSimple linear LIR analysis results \n")
    cat("\nCall:\n")
    print(x.s.linlir$call)
    cat("\nRanges of parameter values of the undominated functions:\n")
    cat("intercept of f in [", eval(x.s.linlir$a.undom[1]), ",", 
        eval(x.s.linlir$a.undom[2]), "]\n", sep = "")
    cat("slope of f in [", eval(x.s.linlir$b.undom[1]), ",", 
        eval(x.s.linlir$b.undom[2]), "]\n", sep = "")
    cat("\nBandwidth: ", eval(2 * x.s.linlir$q.lrm), "\n")
    cat("\nEstimated parameters of the function f.lrm:\n")
    if (!is.vector(x.s.linlir$f.lrm)) {
        cat("intercepts of f.lrm: ", eval(x.s.linlir$f.lrm[, 
            1]), "\n")
        cat("slopes of f.lrm: ", eval(x.s.linlir$f.lrm[, 2]), 
            "\n")
    }
    else {
        cat("intercept of f.lrm: ", eval(x.s.linlir$f.lrm[1]), 
            "\n")
        cat("slope of f.lrm: ", eval(x.s.linlir$f.lrm[2]), "\n")
    }
    cat("\nNumber of observations:", x.s.linlir$n, "\n")
    cat("\nLIR settings:\n")
    cat(paste("p:", eval(x.s.linlir$config$p), "   beta:", eval(x.s.linlir$config$bet), 
        "   epsilon:", eval(x.s.linlir$config$epsilon), "   k.l:", 
        eval(x.s.linlir$config$k.l), "   k.u:", eval(x.s.linlir$config$k.u), 
        " \n"))
    cat(paste("confidence level of each confidence interval:", 
        eval(round(x.s.linlir$config$conf.lev.ci * 100, 2)), 
        "% \n"))
}
undom.a <-
function (dat, b, q.lrm, p = 0.5, bet, epsilon = 0) 
{
    n <- nrow(dat)
    k.l <- kl.ku(n, p, bet, epsilon)[1]
    dat.b <- dat[rowSums(abs(dat)) < Inf, ]
    n.b <- nrow(dat.b)
    dat.b0 <- dat[rowSums(abs(dat[, c(3, 4)])) < Inf, ]
    n.b0 <- nrow(dat.b0)
    d <- matrix(0, n, 2)
    if (b > 0) {
        d[, 1] <- dat[, 3] - b * dat[, 2]
        d[, 2] <- dat[, 4] - b * dat[, 1]
    }
    else {
        if (b == 0) {
            d[, 1] <- dat[, 3]
            d[, 2] <- dat[, 4]
        }
        else {
            d[, 1] <- dat[, 3] - b * dat[, 1]
            d[, 2] <- dat[, 4] - b * dat[, 2]
        }
    }
    d.l <- sort(d[, 1])
    d.u <- sort(d[, 2])
    a.undom <- matrix(0, (n - k.l), 2)
    for (j in 1:nrow(a.undom)) {
        a.undom[j, ] <- c(d.l[k.l + j] - q.lrm, d.u[j] + q.lrm)
    }
    result1 <- a.undom[order(a.undom[, 1]), ]
    attr(result1, "dimnames")[[2]] <- c("a.l", "a.u")
    if (nrow(matrix(a.undom[a.undom[, 1] <= a.undom[, 2], ], 
        ncol = 2)) < 1) {
        result2 <- "There is no undominated line with the chosen slope b!"
    }
    else {
        a.undom <- unique(a.undom, MARGIN = 1)
        preresult2 <- a.undom[a.undom[, 1] <= a.undom[, 2], ]
        if (is.vector(preresult2) == T) {
            result2 <- matrix(preresult2, 1, 2)
            attr(result2, "dimnames")[[2]] <- c("a.l", "a.u")
        }
        else {
            result2 <- preresult2
            i <- 1
            while (i < nrow(result2)) {
                for (k in 1:(nrow(result2) - i)) {
                  if (result2[i, 2] >= result2[i + k, 1]) {
                    result2[i, 2] <- max(result2[i, 2], result2[i + 
                      k, 2])
                  }
                }
                if (result2[i, 2] == max(result2[, 2])) {
                  i <- nrow(result2)
                }
                else {
                  i <- i + 1
                }
            }
            for (i in 2:nrow(result2)) {
                if (result2[i, 2] <= result2[i - 1, 2]) {
                  result2[i, 1] <- result2[i - 1, 1]
                  result2[i, 2] <- result2[i - 1, 2]
                }
            }
            result2 <- unique(result2)
            attr(result2, "dimnames")[[2]] <- c("a.l", "a.u")
        }
    }
    list(result1, result2)
}
undom.para <-
function (dat, b.range, a.grid = 100, q.lrm, p = 0.5, bet, epsilon = 0) 
{
    n <- nrow(dat)
    b.pot <- b.range
    para.undom.inf <- matrix(NA, 1, 3)
    attr(para.undom.inf, "dimnames")[[2]] <- c("a.l", "a.u", 
        "b")
    para.undom <- matrix(NA, 1, 2)
    attr(para.undom, "dimnames")[[2]] <- c("a", "b")
    for (i in 1:length(b.pot)) {
        a.undom <- undom.a(dat, b = b.pot[i], q.lrm, p, bet, 
            epsilon)[[2]]
        if (is.matrix(a.undom)) {
            para.undom.inf <- rbind(para.undom.inf, cbind(a.undom, 
                rep(b.pot[i], times = nrow(a.undom))))
            for (k in 1:nrow(a.undom)) {
                if (round(a.undom[[k, 2]], 10) == round(a.undom[[k, 
                  1]], 10)) {
                  para.undom <- rbind(para.undom, matrix(c(round(a.undom[k, 
                    1], 10), b.pot[i]), nrow = 1, ncol = 2))
                }
                else {
                  a.undom[which(a.undom[, 1] == -Inf), 1] <- -1e+09
                  a.undom[which(a.undom[, 1] == Inf), 1] <- 1e+09
                  a.undom[which(a.undom[, 2] == -Inf), 2] <- -1e+09
                  a.undom[which(a.undom[, 2] == Inf), 2] <- 1e+09
                  if ((a.undom[[k, 2]] - a.undom[[k, 1]]) < (a.grid/20)) {
                    k.grid <- round(a.grid/4)
                    para.undom <- rbind(para.undom, matrix(c(seq(a.undom[[k, 
                      1]], a.undom[[k, 2]], (a.undom[[k, 2]] - 
                      a.undom[[k, 1]])/(k.grid - 1)), rep(b.pot[i], 
                      k.grid)), nrow = k.grid, ncol = 2))
                  }
                  else {
                    if ((a.undom[[k, 2]] - a.undom[[k, 1]]) < 
                      (a.grid/10)) {
                      k.grid <- round(a.grid/2)
                      para.undom <- rbind(para.undom, matrix(c(seq(a.undom[[k, 
                        1]], a.undom[[k, 2]], (a.undom[[k, 2]] - 
                        a.undom[[k, 1]])/(k.grid - 1)), rep(b.pot[i], 
                        k.grid)), nrow = k.grid, ncol = 2))
                    }
                    else {
                      para.undom <- rbind(para.undom, matrix(c(seq(a.undom[[k, 
                        1]], a.undom[[k, 2]], (a.undom[[k, 2]] - 
                        a.undom[[k, 1]])/(a.grid - 1)), rep(b.pot[i], 
                        a.grid)), nrow = a.grid, ncol = 2))
                    }
                  }
                }
            }
        }
    }
    if (which(is.na(para.undom[, 1])) > 1) {
        stop("Too many NA's !\n")
    }
    preresult1 <- para.undom.inf[!is.na(para.undom.inf[, 1]), 
        ]
    preresult2 <- para.undom[!is.na(para.undom[, 1]), ]
    if (!is.vector(preresult1)) {
        result1 <- c(min(preresult1[, 1]), max(preresult1[, 2]))
        attr(result1, "names") <- c("a.min", "a.max")
        if (min(preresult1[, 3]) <= -1e+09 & max(preresult1[preresult1[, 
            3] <= min(preresult1[, 3]), 2]) > 1e+06) {
            result1[2] <- Inf
        }
        if (min(preresult1[, 3]) <= -1e+09 & min(preresult1[preresult1[, 
            3] <= min(preresult1[, 3]), 1]) < -1e+06) {
            result1[1] <- -Inf
        }
        if (max(preresult1[, 3]) >= 1e+09 & min(preresult1[preresult1[, 
            3] >= max(preresult1[, 3]), 1]) < -1e+06) {
            result1[1] <- -Inf
        }
        if (max(preresult1[, 3]) >= 1e+09 & max(preresult1[preresult1[, 
            3] >= max(preresult1[, 3]), 2]) > 1e+06) {
            result1[2] <- Inf
        }
        result2 <- c(min(preresult1[, 3]), max(preresult1[, 3]))
        attr(result2, "names") <- c("b.min", "b.max")
        if (result2[1] <= -1e+09) {
            result2[1] <- -Inf
        }
        if (result2[2] >= 1e+09) {
            result2[2] <- Inf
        }
        result3 <- preresult2
    }
    else {
        result1 <- c(preresult1[1], preresult1[2])
        attr(result1, "names") <- c("a.min", "a.max")
        result2 <- c(preresult1[3], preresult1[3])
        attr(result2, "names") <- c("b.min", "b.max")
        result3 <- preresult2
    }
    list(a.undom = result1, b.undom = result2, undom.para = result3)
}
