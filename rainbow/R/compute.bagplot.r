compute.bagplot = function (x, y, factor = 3, na.rm = FALSE, approx.limit = 300, 
    dkmethod = 2, precision = 1, verbose = FALSE, debug.plots = "no") 
{
    win <- function(dx, dy) {
        atan2(y = dy, x = dx)
    }
    out.of.polygon <- function(xy, pg) {
        if (nrow(pg) == 1) 
            return(xy[, 1] == pg[1] & xy[, 2] == pg[2])
        m <- nrow(xy <- matrix(xy, ncol = 2))
        n <- nrow(pg)
        limit <- -abs(1e-10 * diff(range(pg)))
        pgn <- cbind(diff(c(pg[, 2], pg[1, 2])), -diff(c(pg[, 
            1], pg[1, 1])))
        S <- matrix(colMeans(xy), 1, 2)
        dxy <- cbind(S[1] - pg[, 1], S[2] - pg[, 2])
        S.in.pg <- all(limit < apply(dxy * pgn, 1, sum))
        if (!all(limit < apply(dxy * pgn, 1, sum))) {
            pg <- pg[n:1, ]
            pgn <- -pgn[n:1, ]
        }
        in.pg <- rep(TRUE, m)
        for (j in 1:n) {
            dxy <- xy - matrix(pg[j, ], m, 2, byrow = TRUE)
            in.pg <- in.pg & limit < (dxy %*% pgn[j, ])
        }
        return(!in.pg)
    }
    cut.z.pg <- function(zx, zy, p1x, p1y, p2x, p2y) {
        a2 <- (p2y - p1y) / (p2x - p1x)
        a1 <- zy/zx
        sx <- (p1y - a2 * p1x) / (a1 - a2)
        sy <- a1 * sx
        sxy <- cbind(sx, sy)
        h <- any(is.nan(sxy)) || any(is.na(sxy)) || any(Inf == 
            abs(sxy))
        if (h) {
            if (!exists("verbose")) 
                verbose <- FALSE
            if (verbose) 
                cat("special")
            h <- 0 == (a1 - a2) & sign(zx) == sign(p1x)
            sx <- ifelse(h, p1x, sx)
            sy <- ifelse(h, p1y, sy)
            h <- 0 == (a1 - a2) & sign(zx) != sign(p1x)
            sx <- ifelse(h, p2x, sx)
            sy <- ifelse(h, p2y, sy)
            h <- p1x == p2x & zx != p1x & p1x != 0
            sx <- ifelse(h, p1x, sx)
            sy <- ifelse(h, zy * p1x / zx, sy)
            h <- p1x == p2x & zx != p1x & p1x == 0
            sx <- ifelse(h, p1x, sx)
            sy <- ifelse(h, 0, sy)
            h <- p1x == p2x & zx == p1x & p1x == 0 & sign(zy) == 
                 sign(p1y)
            sx <- ifelse(h, p1x, sx)
            sy <- ifelse(h, p1y, sy)
            h <- p1x == p2x & zx == p1x & p1x == 0 & sign(zy) != 
                 sign(p1y)
            sx <- ifelse(h, p1x, sx)
            sy <- ifelse(h, p2y, sy)
            h <- zx == p1x & zy == p1y
            sx <- ifelse(h, p1x, sx)
            sy <- ifelse(h, p1y, sy)
            h <- zx == p2x & zy == p2y
            sx <- ifelse(h, p2x, sx)
            sy <- ifelse(h, p2y, sy)
            h <- zx == 0 & zy == 0
            sx <- ifelse(h, 0, sx)
            sy <- ifelse(h, 0, sy)
            sxy <- cbind(sx, sy)
        }
        if (!exists("debug.plots")) 
            debug.plots <- "no"
        if (debug.plots == "all") {
            segments(sxy[, 1], sxy[, 2], zx, zy, col = "red")
            segments(0, 0, sxy[, 1], sxy[, 2], type = "l", col = "green", 
                     lty = 2)
            points(sxy, col = "red")
        }
        return(sxy)
    }
    find.cut.z.pg <- function(z, pg, center = c(0, 0), debug.plots = "no") {
        if (!is.matrix(z)) 
            z <- rbind(z)
        if (1 == nrow(pg)) 
            return(matrix(center, nrow(z), 2, TRUE))
        n.pg <- nrow(pg)
        n.z <- nrow(z)
        z <- cbind(z[, 1] - center[1], z[, 2] - center[2])
        pgo <- pg
        pg <- cbind(pg[, 1] - center[1], pg[, 2] - center[2])
        if (!exists("debug.plots")) 
            debug.plots <- "no"
        if (debug.plots == "all") {
            plot(rbind(z, pg, 0), bty = "n")
            points(z, pch = "p")
            lines(c(pg[, 1], pg[1, 1]), c(pg[, 2], pg[1, 2]))
        }
        apg <- win(pg[, 1], pg[, 2])
        apg[is.nan(apg)] <- 0
        a <- order(apg)
        apg <- apg[a]
        pg <- pg[a, ]
        az <- win(z[, 1], z[, 2])
        segm.no <- apply((outer(apg, az, "<")), 2, sum)
        segm.no <- ifelse(segm.no == 0, n.pg, segm.no)
        next.no <- 1 + (segm.no%%length(apg))
        cuts <- cut.z.pg(z[, 1], z[, 2], pg[segm.no, 1], pg[segm.no, 
            2], pg[next.no, 1], pg[next.no, 2])
        cuts <- cbind(cuts[, 1] + center[1], cuts[, 2] + center[2])
        return(cuts)
    }
    hdepth.of.points <- function(tp, n) {
        n.tp <- nrow(tp)
        tphdepth <- rep(0, n.tp)
        dpi <- 2 * pi - 1e-06
        minusplus <- c(rep(-1, n), rep(1, n))
        for (j in 1:n.tp) {
            dx <- tp[j, 1] - xy[, 1]
            dy <- tp[j, 2] - xy[, 2]
            a <- win(dx, dy) + pi
            h <- a < 10
            a <- a[h]
            ident <- sum(!h)
            init <- sum(a < pi)
            a.shift <- (a + pi) %% dpi
            minusplus <- c(rep(-1, length(a)), rep(1, length(a)))
            h <- cumsum(minusplus[order(c(a, a.shift))])
            tphdepth[j] <- init + min(h) + 1
        }
        tphdepth
    }
    expand.hull <- function(pg, k) {
        resolution <- floor(20 * precision)
        pg0 <- xy[hdepth == 1, ]
        pg0 <- pg0[chull(pg0[, 1], pg0[, 2]), ]
        end.points <- find.cut.z.pg(pg, pg0, center = center, 
            debug.plots = debug.plots)
        lam <- ((0 : resolution)^1) / resolution^1
        pg.new <- pg
        for (i in 1:nrow(pg)) {
            tp <- cbind(pg[i, 1] + lam * (end.points[i, 1] - 
                pg[i, 1]), pg[i, 2] + lam * (end.points[i, 2] - 
                pg[i, 2]))
            hd.tp <- hdepth.of.points(tp, nrow(xy))
            ind <- max(sum(hd.tp >= k), 1)
            if (ind < length(hd.tp)) {
                tp <- cbind(tp[ind, 1] + lam * (tp[ind + 1, 1] - 
                  tp[ind, 1]), tp[ind, 2] + lam * (tp[ind + 1, 
                  2] - tp[ind, 2]))
                hd.tp <- hdepth.of.points(tp, nrow(xy))
                ind <- max(sum(hd.tp >= k), 1)
            }
            pg.new[i, ] <- tp[ind, ]
        }
        pg.new <- pg.new[chull(pg.new[, 1], pg.new[, 2]), ]
        pg.add <- 0.5 * (pg.new + rbind(pg.new[-1, ], pg.new[1, 
            ]))
        end.points <- find.cut.z.pg(pg.add, pg0, center = center)
        for (i in 1:nrow(pg.add)) {
            tp <- cbind(pg.add[i, 1] + lam * (end.points[i, 1] - 
                pg.add[i, 1]), pg.add[i, 2] + lam * (end.points[i, 
                2] - pg.add[i, 2]))
            hd.tp <- hdepth.of.points(tp, nrow(xy))
            ind <- max(sum(hd.tp >= k), 1)
            if (ind < length(hd.tp)) {
                tp <- cbind(tp[ind, 1] + lam * (tp[ind + 1, 1] - 
                  tp[ind, 1]), tp[ind, 2] + lam * (tp[ind + 1, 
                  2] - tp[ind, 2]))
                hd.tp <- hdepth.of.points(tp, nrow(xy))
                ind <- max(sum(hd.tp >= k), 1)
            }
            pg.add[i, ] <- tp[ind, ]
        }
        pg.new <- rbind(pg.new, pg.add)
        pg.new <- pg.new[chull(pg.new[, 1], pg.new[, 2]), ]
    }
    cut.p.sl.p.sl <- function(xy1, m1, xy2, m2) {
        sx <- (xy2[2] - m2 * xy2[1] - xy1[2] + m1 * xy1[1])/(m1 - 
               m2)
        sy <- xy1[2] - m1 * xy1[1] + m1 * sx
        if (!is.nan(sy)) 
            return(c(sx, sy))
        if (abs(m1) == Inf) 
            return(c(xy1[1], xy2[2] + m2 * (xy1[1] - xy2[1])))
        if (abs(m2) == Inf) 
            return(c(xy2[1], xy1[2] + m1 * (xy2[1] - xy1[1])))
    }
    pos.to.pg <- function(z, pg, reverse = FALSE) {
        if (reverse) {
            int.no <- apply(outer(pg[, 1], z[, 1], ">="), 2, 
                            sum)
            zy.on.pg <- pg[int.no, 2] + pg[int.no, 3] * (z[, 
                        1] - pg[int.no, 1])
        }
        else {
            int.no <- apply(outer(pg[, 1], z[, 1], "<="), 2, 
                      sum)
            zy.on.pg <- pg[int.no, 2] + pg[int.no, 3] * (z[, 
                        1] - pg[int.no, 1])
        }
        ifelse(z[, 2] < zy.on.pg, "lower", "higher")
    }
    find.polygon.center <- function(xy) {
        if (length(xy) == 2) 
            return(xy[1:2])
        n <- length(xy[, 1])
        mxy <- colMeans(xy)
        xy2 <- rbind(xy[-1, ], xy[1, ])
        xy3 <- cbind(rep(mxy[1], n), mxy[2])
        S <- (xy + xy2 + xy3)/3
        F2 <- abs((xy[, 1] - xy3[, 1]) * (xy2[, 2] - xy3[, 2]) - 
            (xy[, 2] - xy3[, 2]) * (xy2[, 1] - xy3[, 1]))
        lambda <- F2/sum(F2)
        SP <- colSums(cbind(S[, 1] * lambda, S[, 2] * lambda))
        return(SP)
    }
    xydata <- if (missing(y)) 
        x
    else cbind(x, y)
    if (is.data.frame(xydata)) 
        xydata <- as.matrix(xydata)
    if (any(is.na(xydata))) {
        if (na.rm) {
            xydata <- xydata[!apply(is.na(xydata), 1, any), , 
                            drop = FALSE]
            print("Warning: NA elements have been removed!!")
        }
        else {
            xy.means <- colMeans(xydata, na.rm = TRUE)
            for (j in 1:ncol(xydata)) xydata[is.na(xydata[, j]), 
                j] <- xy.means[j]
            print("Warning: NA elements have been exchanged by mean values!!")
        }
    }
    if (nrow(xydata) < 3) {
        print("not enough data points")
        return()
    }
    very.large.data.set <- nrow(xydata) > approx.limit
    set.seed(random.seed <- 13)
    if (very.large.data.set) {
        ind <- sample(seq(nrow(xydata)), size = approx.limit)
        xy <- xydata[ind, ]
    }
    else xy <- xydata
    n <- nrow(xy)
    points.in.bag <- floor(n/2)
    if (verbose) 
        cat("end of initialization")
    prdata <- prcomp(xydata)
    is.one.dim <- (min(prdata[[1]])/max(prdata[[1]])) < 1e-04
    if (is.one.dim) {
        if (verbose) 
            cat("data set one dimensional")
        center <- colMeans(xydata)
        res <- list(xy = xy, xydata = xydata, prdata = prdata, 
            is.one.dim = is.one.dim, center = center)
        class(res) <- "bagplot"
        return(res)
    }
    if (verbose) 
        cat("data not linear")
    xym <- apply(xy, 2, mean)
    xysd <- apply(xy, 2, sd)
    xyxy <- cbind((xy[, 1] - xym[1])/xysd[1], (xy[, 2] - xym[2])/xysd[2])
    dx <- (outer(xy[, 1], xy[, 1], "-"))
    dy <- (outer(xy[, 2], xy[, 2], "-"))
    alpha <- atan2(y = dy, x = dx)
    diag(alpha) <- 1000
    for (j in 1:n) alpha[, j] <- sort(alpha[, j])
    alpha <- alpha[-n, ]
    m <- n - 1
    if (debug.plots == "all") {
        plot(xy, bty = "n")
        xdelta <- abs(diff(range(xy[, 1])))
        dx <- xdelta * 0.3
        for (j in 1:n) {
            p <- xy[j, ]
            dy <- dx * tan(alpha[, j])
            segments(p[1] - dx, p[2] - dy, p[1] + dx, p[2] + 
                dy, col = j)
            text(p[1] - xdelta * 0.02, p[2], j, col = j)
        }
    }
    if (verbose) 
        print("end of computation of angles")
    hdepth <- rep(0, n)
    dpi <- 2 * pi - 1e-06
    mypi <- pi - 1e-06
    minusplus <- c(rep(-1, m), rep(1, m))
    for (j in 1:n) {
        a <- alpha[, j] + pi
        h <- a < 10
        a <- a[h]
        init <- sum(a < mypi)
        a.shift <- (a + pi) %% dpi
        minusplus <- c(rep(-1, length(a)), rep(1, length(a)))
        h <- cumsum(minusplus[order(c(a, a.shift))])
        hdepth[j] <- init + min(h) + 1
    }
    if (verbose) {
        print("end of computation of hdepth:")
        print(hdepth)
    }
    if (debug.plots == "all") {
        plot(xy, bty = "n")
        xdelta <- abs(diff(range(xy[, 1])))
        dx <- xdelta * 0.1
        for (j in 1:n) {
            a <- alpha[, j] + pi
            a <- a[a < 10]
            init <- sum(a < pi)
            a.shift <- (a + pi) %% dpi
            minusplus <- c(rep(-1, length(a)), rep(1, length(a)))
            h <- cumsum(minusplus[ao <- (order(c(a, a.shift)))])
            no <- which((init + min(h)) == (init + h))[1]
            p <- xy[j, ]
            dy <- dx * tan(alpha[, j])
            segments(p[1] - dx, p[2] - dy, p[1] + dx, p[2] + 
                     dy, col = j, lty = 3)
            dy <- dx * tan(c(sort(a), sort(a))[no])
            segments(p[1] - 5 * dx, p[2] - 5 * dy, p[1] + 5 * 
                     dx, p[2] + 5 * dy, col = "black")
            text(p[1] - xdelta * 0.02, p[2], hdepth[j], col = 1, 
                     cex = 2.5)
        }
    }
    hd.table <- table(sort(hdepth))
    d.k <- cbind(dk = rev(cumsum(rev(hd.table))), k = as.numeric(names(hd.table)))
    k.1 <- sum(points.in.bag < d.k[, 1])
    if (nrow(d.k) > k.1) {
        k <- d.k[k.1 + 1, 2]
    }
    else {
        k <- d.k[k.1, 2]
    }
    if (verbose) {
        cat("numbers of members of dk:")
        print(hd.table)
    }
    if (verbose) {
        cat("end of computation of k, k=", k)
    }
    center <- apply(xy[which(hdepth == max(hdepth)), , drop = FALSE], 
                    2, mean)
    hull.center <- NULL
    if (5 < nrow(xy) && length(hd.table) > 2) {
        n.p <- floor(c(32, 16, 8)[1 + (n > 50) + (n > 200)] * 
                     precision)
        h <- cands <- xy[rev(order(hdepth))[1:6], ]
        cands <- cands[chull(cands[, 1], cands[, 2]), ]
        n.c <- nrow(cands)
        if (is.null(n.c)) 
            cands <- h
        xyextr <- rbind(apply(cands, 2, min), apply(cands, 2, 
                        max))
        xydel <- 2 * (xyextr[2, ] - xyextr[1, ])/n.p
        h1 <- seq(xyextr[1, 1], xyextr[2, 1], length = n.p)
        h2 <- seq(xyextr[1, 2], xyextr[2, 2], length = n.p)
        tp <- cbind(matrix(h1, n.p, n.p)[1:n.p^2], matrix(h2, 
                    n.p, n.p, TRUE)[1:n.p^2])
        tphdepth <- max(hdepth.of.points(tp, n)) - 1
        num <- floor(c(417, 351, 171, 85, 67, 43)[sum(n > c(1, 
                     50, 100, 150, 200, 250))] * precision)
        num.h <- floor(num/2)
        angles <- seq(0, pi, length = num.h)
        ang <- tan(pi/2 - angles)
        kkk <- tphdepth
        if (verbose) {
            cat("max-hdepth found:")
            print(kkk)
        }
        ia <- 1
        a <- angles[ia]
        xyt <- xyxy %*% c(cos(a), -sin(a))
        xyto <- order(xyt)
        ind.k <- xyto[kkk]
        cutp <- c(xyxy[ind.k, 1], -10)
        dxy <- diff(range(xyxy))
        pg <- rbind(c(cutp[1], -dxy, Inf), c(cutp[1], dxy, NA))
        ind.kk <- xyto[n + 1 - kkk]
        cutpl <- c(xyxy[ind.kk, 1], 10)
        pgl <- rbind(c(cutpl[1], dxy, Inf), c(cutpl[1], -dxy, 
            NA))
        if (debug.plots == "all") {
            plot(xyxy, type = "p", bty = "n")
        }
        for (ia in seq(angles)[-1]) {
            a <- angles[ia]
            angtan <- ang[ia]
            xyt <- xyxy %*% c(cos(a), -sin(a))
            xyto <- order(xyt)
            ind.k <- xyto[kkk]
            ind.kk <- xyto[n + 1 - kkk]
            pnew <- xyxy[ind.k, ]
            pnewl <- xyxy[ind.kk, ]
            if (debug.plots == "all") 
                points(pnew[1], pnew[2], col = "red")
            if (abs(angtan) > 1e+10) {
                pg.no <- sum(pg[, 1] < pnew[1])
                cutp <- c(pnew[1], pg[pg.no, 2] + pg[pg.no, 3] * 
                        (pnew[1] - pg[pg.no, 1]))
                pg.nol <- sum(pgl[, 1] >= pnewl[1])
                cutpl <- c(pnewl[1], pgl[pg.nol, 2] + pgl[pg.nol, 
                         3] * (pnewl[1] - pgl[pg.nol, 1]))
            }
            else {
                pg.inter <- pg[, 2] - angtan * pg[, 1]
                pnew.inter <- pnew[2] - angtan * pnew[1]
                pg.no <- sum(pg.inter < pnew.inter)
                cutp <- cut.p.sl.p.sl(pnew, ang[ia], pg[pg.no, 
                         1:2], pg[pg.no, 3])
                pg.interl <- pgl[, 2] - angtan * pgl[, 1]
                pnew.interl <- pnewl[2] - angtan * pnewl[1]
                pg.nol <- sum(pg.interl > pnew.interl)
                cutpl <- cut.p.sl.p.sl(pnewl, angtan, pgl[pg.nol, 
                         1:2], pgl[pg.nol, 3])
            }
            pg <- rbind(pg[1:pg.no, ], c(cutp, angtan), c(cutp[1] + 
                dxy, cutp[2] + angtan * dxy, NA))
            pgl <- rbind(pgl[1:pg.nol, ], c(cutpl, angtan), c(cutpl[1] - 
                dxy, cutpl[2] - angtan * dxy, NA))
            if (debug.plots == "all") {
                points(pnew[1], pnew[2], col = "red")
                hx <- xyxy[ind.k, c(1, 1)]
                hy <- xyxy[ind.k, c(2, 2)]
                segments(hx, hy, c(10, -10), hy + ang[ia] * (c(10, 
                  -10) - hx), lty = 2)
                points(cutpl[1], cutpl[2], col = "red")
                hx <- xyxy[ind.kk, c(1, 1)]
                hy <- xyxy[ind.kk, c(2, 2)]
                segments(hx, hy, c(10, -10), hy + ang[ia] * (c(10, 
                  -10) - hx), lty = 2)
            }
        }
        limit <- 1e-10
        pg <- pg[c(TRUE, (abs(diff(pg[, 1])) > limit) | (abs(diff(pg[, 
            2])) > limit)), ]
        pgl <- pgl[c(TRUE, (abs(diff(pgl[, 1])) > limit) | (abs(diff(pgl[, 
            2])) > limit)), ]
        pg <- pg[-nrow(pg), ][-1, , drop = FALSE]
        pgl <- pgl[-nrow(pgl), ][-1, , drop = FALSE]
        indl <- pos.to.pg(pgl, pg)
        indu <- pos.to.pg(pg, pgl, TRUE)
        sr <- sl <- NULL
        if (indu[(npg <- nrow(pg))] == "lower" & indl[1] == "higher") {
            rnuml <- which(indl == "lower")[1] - 1
            rnumu <- npg + 1 - which(rev(indu == "higher"))[1]
            if (is.na(rnuml)) 
                rnuml <- sum(pg[rnumu, 1] < pgl[, 1])
            if (is.na(rnumu)) 
                rnumu <- sum(pg[, 1] < pgl[rnuml, 1])
            xyl <- pgl[rnuml, ]
            xyu <- pg[rnumu, ]
            sr <- cut.p.sl.p.sl(xyl[1:2], xyl[3], xyu[1:2], xyu[3])
        }
        if (indl[(npgl <- nrow(pgl))] == "higher" & indu[1] == 
            "lower") {
            lnuml <- npgl + 1 - which(rev(indl == "lower"))[1]
            lnumu <- which(indu == "higher")[1] - 1
            if (is.na(lnuml)) 
                lnuml <- sum(pg[lnumu, 1] < pgl[, 1])
            if (is.na(lnumu)) 
                lnumu <- sum(pg[, 1] < pgl[lnuml, 1])
            xyl <- pgl[lnuml, ]
            xyu <- pg[lnumu, ]
            sl <- cut.p.sl.p.sl(xyl[1:2], xyl[3], xyu[1:2], xyu[3])
        }
        pg <- rbind(pg[indu == "higher", 1:2, drop = FALSE], 
                    sr, pgl[indl == "lower", 1:2, drop = FALSE], sl)
        if (debug.plots == "all") 
            lines(rbind(pg, pg[1, ]), col = "red")
        pg <- pg[chull(pg[, 1], pg[, 2]), ]
        hull.center <- cbind(pg[, 1] * xysd[1] + xym[1], pg[, 
            2] * xysd[2] + xym[2])
        center <- find.polygon.center(hull.center)
        if (verbose) {
            cat("hull.center", hull.center)
            print(table(tphdepth))
        }
    }
    if (verbose) 
        cat("center depth:", hdepth.of.points(rbind(center), 
            n) - 1)
    if (verbose) {
        print("end of computation of center")
        print(center)
    }
    if (dkmethod == 1) {
        xyi <- xy[hdepth >= k, , drop = FALSE]
        pdk <- xyi[chull(xyi[, 1], xyi[, 2]), , drop = FALSE]
        xyo <- xy[hdepth >= (k - 1), , drop = FALSE]
        pdk.1 <- xyo[chull(xyo[, 1], xyo[, 2]), , drop = FALSE]
        if (verbose) 
            cat("hull computed:")
        if (debug.plots == "all") {
            plot(xy, bty = "n")
            h <- rbind(pdk, pdk[1, ])
            lines(h, col = "red", lty = 2)
            h <- rbind(pdk.1, pdk.1[1, ])
            lines(h, col = "blue", lty = 3)
            points(center[1], center[2], pch = 8, col = "red")
        }
        exp.dk <- expand.hull(pdk, k)
        exp.dk.1 <- expand.hull(exp.dk, k - 1)
    }
    else {
        num <- floor(c(417, 351, 171, 85, 67, 43)[sum(n > c(1, 
            50, 100, 150, 200, 250))] * precision)
        num.h <- floor(num/2)
        angles <- seq(0, pi, length = num.h)
        ang <- tan(pi/2 - angles)
        kkk <- k
        ia <- 1
        a <- angles[ia]
        xyt <- xyxy %*% c(cos(a), -sin(a))
        xyto <- order(xyt)
        ind.k <- xyto[kkk]
        cutp <- c(xyxy[ind.k, 1], -10)
        dxy <- diff(range(xyxy))
        pg <- rbind(c(cutp[1], -dxy, Inf), c(cutp[1], dxy, NA))
        ind.kk <- xyto[n + 1 - kkk]
        cutpl <- c(xyxy[ind.kk, 1], 10)
        pgl <- rbind(c(cutpl[1], dxy, Inf), c(cutpl[1], -dxy, 
            NA))
        if (debug.plots == "all") {
            plot(xyxy, type = "p", bty = "n")
        }
        for (ia in seq(angles)[-1]) {
            a <- angles[ia]
            angtan <- ang[ia]
            xyt <- xyxy %*% c(cos(a), -sin(a))
            xyto <- order(xyt)
            ind.k <- xyto[kkk]
            ind.kk <- xyto[n + 1 - kkk]
            pnew <- xyxy[ind.k, ]
            pnewl <- xyxy[ind.kk, ]
            if (debug.plots == "all") 
                points(pnew[1], pnew[2], col = "red")
            if (abs(angtan) > 1e+10) {
                pg.no <- sum(pg[, 1] < pnew[1])
                cutp <- c(pnew[1], pg[pg.no, 2] + pg[pg.no, 3] * 
                  (pnew[1] - pg[pg.no, 1]))
                pg.nol <- sum(pgl[, 1] >= pnewl[1])
                cutpl <- c(pnewl[1], pgl[pg.nol, 2] + pgl[pg.nol, 
                  3] * (pnewl[1] - pgl[pg.nol, 1]))
            }
            else {
                pg.inter <- pg[, 2] - angtan * pg[, 1]
                pnew.inter <- pnew[2] - angtan * pnew[1]
                pg.no <- sum(pg.inter < pnew.inter)
                cutp <- cut.p.sl.p.sl(pnew, ang[ia], pg[pg.no, 
                  1:2], pg[pg.no, 3])
                pg.interl <- pgl[, 2] - angtan * pgl[, 1]
                pnew.interl <- pnewl[2] - angtan * pnewl[1]
                pg.nol <- sum(pg.interl > pnew.interl)
                cutpl <- cut.p.sl.p.sl(pnewl, angtan, pgl[pg.nol, 
                  1:2], pgl[pg.nol, 3])
            }
            pg <- rbind(pg[1:pg.no, ], c(cutp, angtan), c(cutp[1] + 
                dxy, cutp[2] + angtan * dxy, NA))
            pgl <- rbind(pgl[1:pg.nol, ], c(cutpl, angtan), c(cutpl[1] - 
                dxy, cutpl[2] - angtan * dxy, NA))
            if (debug.plots == "all") {
                points(pnew[1], pnew[2], col = "red")
                hx <- xyxy[ind.k, c(1, 1)]
                hy <- xyxy[ind.k, c(2, 2)]
                segments(hx, hy, c(10, -10), hy + ang[ia] * (c(10, 
                  -10) - hx), lty = 2)
                points(cutpl[1], cutpl[2], col = "red")
                hx <- xyxy[ind.kk, c(1, 1)]
                hy <- xyxy[ind.kk, c(2, 2)]
                segments(hx, hy, c(10, -10), hy + ang[ia] * (c(10, 
                  -10) - hx), lty = 2)
            }
        }
        limit <- 1e-10
        pg <- pg[c(TRUE, (abs(diff(pg[, 1])) > limit) | (abs(diff(pg[, 
            2])) > limit)), ]
        pgl <- pgl[c(TRUE, (abs(diff(pgl[, 1])) > limit) | (abs(diff(pgl[, 
            2])) > limit)), ]
        pg <- pg[-nrow(pg), ][-1, , drop = FALSE]
        pgl <- pgl[-nrow(pgl), ][-1, , drop = FALSE]
        indl <- pos.to.pg(pgl, pg)
        indu <- pos.to.pg(pg, pgl, TRUE)
        sr <- sl <- NULL
        if (indu[(npg <- nrow(pg))] == "lower" & indl[1] == "higher") {
            rnuml <- which(indl == "lower")[1] - 1
            rnumu <- npg + 1 - which(rev(indu == "higher"))[1]
            if (is.na(rnuml)) 
                rnuml <- sum(pg[rnumu, 1] < pgl[, 1])
            if (is.na(rnumu)) 
                rnumu <- sum(pg[, 1] < pgl[rnuml, 1])
            xyl <- pgl[rnuml, ]
            xyu <- pg[rnumu, ]
            sr <- cut.p.sl.p.sl(xyl[1:2], xyl[3], xyu[1:2], xyu[3])
        }
        if (indl[(npgl <- nrow(pgl))] == "higher" & indu[1] == 
            "lower") {
            lnuml <- npgl + 1 - which(rev(indl == "lower"))[1]
            lnumu <- which(indu == "higher")[1] - 1
            if (is.na(lnuml)) 
                lnuml <- sum(pg[lnumu, 1] < pgl[, 1])
            if (is.na(lnumu)) 
                lnumu <- sum(pg[, 1] < pgl[lnuml, 1])
            xyl <- pgl[lnuml, ]
            xyu <- pg[lnumu, ]
            sl <- cut.p.sl.p.sl(xyl[1:2], xyl[3], xyu[1:2], xyu[3])
        }
        pg <- rbind(pg[indu == "higher", 1:2, drop = FALSE], 
            sr, pgl[indl == "lower", 1:2, drop = FALSE], sl)
        if (debug.plots == "all") 
            lines(rbind(pg, pg[1, ]), col = "red")
        pg <- pg[chull(pg[, 1], pg[, 2]), ]
        exp.dk <- cbind(pg[, 1] * xysd[1] + xym[1], pg[, 2] * 
                        xysd[2] + xym[2])
        if (kkk > 1) 
            kkk <- kkk - 1
        ia <- 1
        a <- angles[ia]
        xyt <- xyxy %*% c(cos(a), -sin(a))
        xyto <- order(xyt)
        ind.k <- xyto[kkk]
        cutp <- c(xyxy[ind.k, 1], -10)
        dxy <- diff(range(xyxy))
        pg <- rbind(c(cutp[1], -dxy, Inf), c(cutp[1], dxy, NA))
        ind.kk <- xyto[n + 1 - kkk]
        cutpl <- c(xyxy[ind.kk, 1], 10)
        pgl <- rbind(c(cutpl[1], dxy, Inf), c(cutpl[1], -dxy, 
            NA))
        if (debug.plots == "all") {
            plot(xyxy, type = "p", bty = "n")
        }
        for (ia in seq(angles)[-1]) {
            a <- angles[ia]
            angtan <- ang[ia]
            xyt <- xyxy %*% c(cos(a), -sin(a))
            xyto <- order(xyt)
            ind.k <- xyto[kkk]
            ind.kk <- xyto[n + 1 - kkk]
            pnew <- xyxy[ind.k, ]
            pnewl <- xyxy[ind.kk, ]
            if (debug.plots == "all") 
                points(pnew[1], pnew[2], col = "red")
            if (abs(angtan) > 1e+10) {
                pg.no <- sum(pg[, 1] < pnew[1])
                cutp <- c(pnew[1], pg[pg.no, 2] + pg[pg.no, 3] * 
                  (pnew[1] - pg[pg.no, 1]))
                pg.nol <- sum(pgl[, 1] >= pnewl[1])
                cutpl <- c(pnewl[1], pgl[pg.nol, 2] + pgl[pg.nol, 
                  3] * (pnewl[1] - pgl[pg.nol, 1]))
            }
            else {
                pg.inter <- pg[, 2] - angtan * pg[, 1]
                pnew.inter <- pnew[2] - angtan * pnew[1]
                pg.no <- sum(pg.inter < pnew.inter)
                cutp <- cut.p.sl.p.sl(pnew, ang[ia], pg[pg.no, 
                  1:2], pg[pg.no, 3])
                pg.interl <- pgl[, 2] - angtan * pgl[, 1]
                pnew.interl <- pnewl[2] - angtan * pnewl[1]
                pg.nol <- sum(pg.interl > pnew.interl)
                cutpl <- cut.p.sl.p.sl(pnewl, angtan, pgl[pg.nol, 
                  1:2], pgl[pg.nol, 3])
            }
            pg <- rbind(pg[1:pg.no, ], c(cutp, angtan), c(cutp[1] + 
                dxy, cutp[2] + angtan * dxy, NA))
            pgl <- rbind(pgl[1:pg.nol, ], c(cutpl, angtan), c(cutpl[1] - 
                dxy, cutpl[2] - angtan * dxy, NA))
            if (debug.plots == "all") {
                points(pnew[1], pnew[2], col = "red")
                hx <- xyxy[ind.k, c(1, 1)]
                hy <- xyxy[ind.k, c(2, 2)]
                segments(hx, hy, c(10, -10), hy + ang[ia] * (c(10, 
                  -10) - hx), lty = 2)
                points(cutpl[1], cutpl[2], col = "red")
                hx <- xyxy[ind.kk, c(1, 1)]
                hy <- xyxy[ind.kk, c(2, 2)]
                segments(hx, hy, c(10, -10), hy + ang[ia] * (c(10, 
                  -10) - hx), lty = 2)
            }
        }
        limit <- 1e-10
        pg <- pg[c(TRUE, (abs(diff(pg[, 1])) > limit) | (abs(diff(pg[, 
            2])) > limit)), ]
        pgl <- pgl[c(TRUE, (abs(diff(pgl[, 1])) > limit) | (abs(diff(pgl[, 
            2])) > limit)), ]
        pg <- pg[-nrow(pg), ][-1, , drop = FALSE]
        pgl <- pgl[-nrow(pgl), ][-1, , drop = FALSE]
        indl <- pos.to.pg(pgl, pg)
        indu <- pos.to.pg(pg, pgl, TRUE)
        sr <- sl <- NULL
        if (indu[(npg <- nrow(pg))] == "lower" & indl[1] == "higher") {
            rnuml <- which(indl == "lower")[1] - 1
            rnumu <- npg + 1 - which(rev(indu == "higher"))[1]
            if (is.na(rnuml)) 
                rnuml <- sum(pg[rnumu, 1] < pgl[, 1])
            if (is.na(rnumu)) 
                rnumu <- sum(pg[, 1] < pgl[rnuml, 1])
            xyl <- pgl[rnuml, ]
            xyu <- pg[rnumu, ]
            sr <- cut.p.sl.p.sl(xyl[1:2], xyl[3], xyu[1:2], xyu[3])
        }
        if (indl[(npgl <- nrow(pgl))] == "higher" & indu[1] == 
            "lower") {
            lnuml <- npgl + 1 - which(rev(indl == "lower"))[1]
            lnumu <- which(indu == "higher")[1] - 1
            if (is.na(lnuml)) 
                lnuml <- sum(pg[lnumu, 1] < pgl[, 1])
            if (is.na(lnumu)) 
                lnumu <- sum(pg[, 1] < pgl[lnuml, 1])
            xyl <- pgl[lnuml, ]
            xyu <- pg[lnumu, ]
            sl <- cut.p.sl.p.sl(xyl[1:2], xyl[3], xyu[1:2], xyu[3])
        }
        pg <- rbind(pg[indu == "higher", 1:2, drop = FALSE], 
            sr, pgl[indl == "lower", 1:2, drop = FALSE], sl)
        if (debug.plots == "all") 
            lines(rbind(pg, pg[1, ]), col = "red")
        pg <- pg[chull(pg[, 1], pg[, 2]), ]
        exp.dk.1 <- cbind(pg[, 1] * xysd[1] + xym[1], pg[, 2] * 
            xysd[2] + xym[2])
    }
    if (nrow(d.k) == k.1 || nrow(d.k) == 1) 
        lambda <- 0
    else {
        lambda <- (n/2 - d.k[k.1 + 1, 1])/(d.k[k.1, 1] - d.k[k.1 + 
            1, 1])
    }
    if (verbose) 
        cat("lambda", lambda)
    cut.on.pdk.1 <- find.cut.z.pg(exp.dk, exp.dk.1, center = center)
    cut.on.pdk <- find.cut.z.pg(exp.dk.1, exp.dk, center = center)
    h1 <- (1 - lambda) * exp.dk + lambda * cut.on.pdk.1
    h2 <- (1 - lambda) * cut.on.pdk + lambda * exp.dk.1
    h <- rbind(h1, h2)
    h <- h[!is.nan(h[, 1]) & !is.nan(h[, 2]), ]
    hull.bag <- h[chull(h[, 1], h[, 2]), ]
    if (verbose) 
        cat("bag completed:")
    if (debug.plots == "all") {
        lines(hull.bag, col = "red")
    }
    hull.loop <- cbind(hull.bag[, 1] - center[1], hull.bag[, 
        2] - center[2])
    hull.loop <- factor * hull.loop
    hull.loop <- cbind(hull.loop[, 1] + center[1], hull.loop[, 
        2] + center[2])
    if (verbose) 
        cat("loop computed")
    if (!very.large.data.set) {
        pxy.bag <- xydata[hdepth >= k, , drop = FALSE]
        pkt.cand <- xydata[hdepth == (k - 1), , drop = FALSE]
        pkt.not.bag <- xydata[hdepth < (k - 1), , drop = FALSE]
        if (length(pkt.cand) > 0) {
            outside <- out.of.polygon(pkt.cand, hull.bag)
            if (sum(!outside) > 0) 
                pxy.bag <- rbind(pxy.bag, pkt.cand[!outside, 
                  ])
            if (sum(outside) > 0) 
                pkt.not.bag <- rbind(pkt.not.bag, pkt.cand[outside, 
                  ])
        }
    }
    else {
        extr <- out.of.polygon(xydata, hull.bag)
        pxy.bag <- xydata[!extr, ]
        pkt.not.bag <- xydata[extr, , drop = FALSE]
    }
    if (length(pkt.not.bag) > 0) {
        extr <- out.of.polygon(pkt.not.bag, hull.loop)
        pxy.outlier <- pkt.not.bag[extr, , drop = FALSE]
        if (0 == length(pxy.outlier)) 
            pxy.outlier <- NULL
        pxy.outer <- pkt.not.bag[!extr, , drop = FALSE]
    }
    else {
        pxy.outer <- pxy.outlier <- NULL
    }
    if (verbose) 
        cat("points of bag, outer points and outlier identified")
    hull.loop <- rbind(pxy.outer, hull.bag)
    hull.loop <- hull.loop[chull(hull.loop[, 1], hull.loop[, 
        2]), ]
    if (verbose) 
        cat("end of computation of loop")
    res <- list(center = center, hull.center = hull.center, hull.bag = hull.bag, 
        hull.loop = hull.loop, pxy.bag = pxy.bag, pxy.outer = if (length(pxy.outer) > 
            0) pxy.outer else NULL, pxy.outlier = if (length(pxy.outlier) > 
            0) pxy.outlier else NULL, hdepths = hdepth, is.one.dim = is.one.dim, 
        prdata = prdata, random.seed = random.seed, xy = xy, 
        xydata = xydata)
    if (verbose) 
        res <- c(res, list(exp.dk = exp.dk, exp.dk.1 = exp.dk.1, 
            hdepth = hdepth))
    class(res) <- "bagplot"
    return(res)
}
