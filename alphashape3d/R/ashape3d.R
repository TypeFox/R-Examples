ashape3d <-
function (x, alpha, pert = FALSE, eps = 1e-09) 
{
    flag <- 1
    flag2 <- 0
    alphaux <- alpha
    if (any(alphaux < 0)) {
        stop("Parameter alpha must be greater or equal to zero", 
            call. = TRUE)
    }
    if (inherits(x, "ashape3d")) {
        inh <- 1
        tc.def <- x$tetra
        tri.def <- x$triang
        ed.def <- x$edge
        vt.def <- x$vertex
        alphaold <- x$alpha
        if (any(match(alphaux, alphaold, nomatch = 0))) {
            warning("Some values of alpha have already been computed", 
                call. = FALSE)
            alphaux <- alphaux[is.na(match(alphaux, alphaold))]
        }
        if (!is.null(x$xpert)) {
            flag <- 1
            xorig <- x$x
            x <- x$xpert
        }
        else {
            x <- x$x
        }
    }
    else {
        if (any(duplicated(x))) {
            warning("Duplicate poits were removed", call. = FALSE, 
                immediate. = TRUE)
            x <- unique(x)
        }
        inh <- 0
        x = x * 1
        xorig <- x
        while (flag == 1) {
            flag <- 0
            n <- dim(x)[1]
            x4 <- x[, 1]^2 + x[, 2]^2 + x[, 3]^2
            tc <- matrix(delaunayn(x), ncol = 4)
            ntc <- dim(tc)[1]
            ORDM <- matrix(as.integer(0), ntc, 4)
            storage.mode(tc) <- "integer"
            storage.mode(ORDM) <- "integer"
            tc <- .C("sortm", ORDM, as.integer(ntc), as.integer(4), 
                tc, PACKAGE = "alphashape3d")[[1]]
            tc.aux <- matrix(1:(4 * ntc), ncol = 4)
            tri <- rbind(tc[, -4], tc[, -3], tc[, -2], tc[, -1])
            t1 <- tri[, 1]
            t2 <- tri[, 2]
            t3 <- tri[, 3]
            ix3 = .Call("sortbycolumn", t1, t2, t3)
            d1 <- abs(diff(t1[ix3]))
            d2 <- abs(diff(t2[ix3]))
            d3 <- abs(diff(t3[ix3]))
            dup <- (d1 + d2 + d3) == 0
            dup <- c(FALSE, dup)
            i1 <- ix3[dup]
            i2 <- ix3[c(dup[-1], FALSE)]
            ih <- (1:(4 * ntc))[-c(i1, i2)]
            ntris <- length(i1)
            ntrich <- length(ih)
            ntri <- ntris + ntrich
            ind.tri <- numeric(4 * ntc)
            ind.tri[i1] <- 1:ntris
            ind.tri[i2] <- 1:ntris
            ind.tri[ih] <- ((ntris + 1):ntri)
            in.tc <- rep(1:ntc, 4)
            t1u <- t1[c(i1, ih)]
            t2u <- t2[c(i1, ih)]
            t3u <- t3[c(i1, ih)]
            on.ch3 <- numeric(ntri)
            on.ch3[(ntris + 1):ntri] <- 1
            tri.aux <- matrix(1:(3 * ntri), ncol = 3)
            storage.mode(x) <- "double"
            storage.mode(t1u) <- "integer"
            storage.mode(t2u) <- "integer"
            storage.mode(t3u) <- "integer"
            m123 <- numeric(ntri)
            storage.mode(m123) <- "double"
            fm123 <- .C("fm123", x, as.integer(n), t1u, t2u, 
                t3u, as.integer(ntri), m123, PACKAGE = "alphashape3d")
            m123 <- fm123[[7]]
            m1230 <- +m123[ind.tri[tc.aux[, 1]]] - m123[ind.tri[tc.aux[, 
                2]]] + m123[ind.tri[tc.aux[, 3]]] - m123[ind.tri[tc.aux[, 
                4]]]
            if (any(m1230 == 0) & !pert) {
                stop("The general position assumption is not satisfied\nPlease enter data in general position or set pert=TRUE to allow perturbation", 
                  call. = FALSE)
            }
            if (any(m1230 == 0) & pert) {
                flag <- 1
            }
            if (flag == 0) {
                e1 <- c(t1u, t1u, t2u)
                e2 <- c(t2u, t3u, t3u)
                a1 <- 10^nchar(max(c(e1, e2)))
                a2 <- e2 + e1 * a1
                ix2 <- order(a2)
                d1 <- abs(diff(e1[ix2]))
                d2 <- abs(diff(e2[ix2]))
                dup <- (d1 + d2) == 0
                in.tri <- rep(1:ntri, 3)
                dup <- c(which(!dup), length(dup) + 1)
                dup.aux <- dup
                dup <- c(dup[1], diff(dup))
                ned <- length(dup)
                auxi <- rep(1:ned, dup)
                ind.ed <- numeric()
                ind.ed[ix2] <- auxi
                in.trio <- in.tri[ix2]
                e1u <- e1[ix2][dup.aux]
                e2u <- e2[ix2][dup.aux]
                on.ch2 <- numeric(ned)
                on.ch2[ind.ed[tri.aux[as.logical(on.ch3), ]]] <- 1
                vt <- c(e1u, e2u)
                ind.vt <- vt
                nvt <- n
                in.ed <- rep(1:ned, 2)
                ix1 <- order(vt)
                dup1 <- diff(vt[ix1]) == 0
                dup1 <- c(which(!dup1), length(dup1) + 1)
                dup1.aux <- dup1
                dup1 <- c(dup1[1], diff(dup1))
                in.edo <- in.ed[ix1]
                on.ch1 <- numeric(nvt)
                on.ch1[c(t1u[as.logical(on.ch3)], t2u[as.logical(on.ch3)], 
                  t3u[as.logical(on.ch3)])] <- 1
                storage.mode(e1u) <- "integer"
                storage.mode(e2u) <- "integer"
                mk0 <- matrix(0, nrow = ned, ncol = 3)
                storage.mode(mk0) <- "double"
                num2 <- numeric(ned)
                storage.mode(num2) <- "double"
                fmk0 <- .C("fmk0", x, as.integer(n), e1u, e2u, 
                  as.integer(ned), mk0, num2, PACKAGE = "alphashape3d")
                mk0 <- fmk0[[6]]
                num.rho2 <- fmk0[[7]]
                rho2 <- sqrt(0.25 * num.rho2)
                storage.mode(tri.aux) <- "integer"
                storage.mode(ind.ed) <- "integer"
                storage.mode(mk0) <- "double"
                m23 <- numeric(ntri)
                storage.mode(m23) <- "double"
                m13 <- numeric(ntri)
                storage.mode(m13) <- "double"
                m12 <- numeric(ntri)
                storage.mode(m12) <- "double"
                storage.mode(num2) <- "double"
                num3 <- numeric(ntri)
                storage.mode(num3) <- "double"
                fmij0 <- .C("fmij0", x, as.integer(n), t1u, t2u, 
                  t3u, as.integer(ntri), tri.aux, ind.ed, as.integer(ned), 
                  mk0, m23, m13, m12, num.rho2, num3, PACKAGE = "alphashape3d")
                m230 <- fmij0[[11]]
                m130 <- fmij0[[12]]
                m120 <- fmij0[[13]]
                num.rho3 <- fmij0[[15]]
                den.rho3 <- m230^2 + m130^2 + m120^2
                if (any(den.rho3 == 0) & !pert) {
                  stop("The general position assumption is not satisfied\nPlease enter data in general position or set pert=TRUE to allow perturbation", 
                    call. = FALSE)
                }
                if (any(den.rho3 == 0) & pert) {
                  flag <- 1
                }
            }
            if (flag == 0) {
                rho3 <- sqrt(0.25 * num.rho3/(m230^2 + m130^2 + 
                  m120^2))
                storage.mode(x4) <- "double"
                storage.mode(tc.aux) <- "integer"
                storage.mode(ind.tri) <- "integer"
                storage.mode(m230) <- "double"
                storage.mode(m130) <- "double"
                storage.mode(m120) <- "double"
                m2340 <- numeric(ntc)
                m1340 <- numeric(ntc)
                m1240 <- numeric(ntc)
                storage.mode(m2340) <- "double"
                storage.mode(m1340) <- "double"
                storage.mode(m1240) <- "double"
                fmijk0 <- .C("fmijk0", x4, as.integer(n), tc, 
                  as.integer(ntc), tc.aux, ind.tri, as.integer(ntri), 
                  m230, m130, m120, m2340, m1340, m1240, PACKAGE = "alphashape3d")
                m2340 <- fmijk0[[11]]
                m1340 <- fmijk0[[12]]
                m1240 <- fmijk0[[13]]
                m1234 <- -(-x4[tc[, 4]] * m123[ind.tri[tc.aux[, 
                  1]]] + x4[tc[, 3]] * m123[ind.tri[tc.aux[, 
                  2]]] - x4[tc[, 2]] * m123[ind.tri[tc.aux[, 
                  3]]] + x4[tc[, 1]] * m123[ind.tri[tc.aux[, 
                  4]]])
                rho.sq <- (0.25 * (m2340^2 + m1340^2 + m1240^2 + 
                  4 * m1230 * m1234)/m1230^2)
                if (any(rho.sq < 0)) {
                  ind <- which(rho.sq < 0)
                  v1 <- x[tc[, 1], ]
                  v2 <- x[tc[, 2], ]
                  v3 <- x[tc[, 3], ]
                  v4 <- x[tc[, 4], ]
                  v1 <- v1[ind, , drop = FALSE]
                  v2 <- v2[ind, , drop = FALSE]
                  v3 <- v3[ind, , drop = FALSE]
                  v4 <- v4[ind, , drop = FALSE]
                  nv1 <- rowSums(v1^2)
                  nv2 <- rowSums(v2^2)
                  nv3 <- rowSums(v3^2)
                  nv4 <- rowSums(v4^2)
                  Dx <- numeric()
                  Dy <- numeric()
                  Dz <- numeric()
                  ct <- numeric()
                  a <- numeric()
                  for (i in 1:length(ind)) {
                    tetra <- rbind(v1[i, ], v2[i, ], v3[i, ], 
                      v4[i, ])
                    nor <- c(nv1[i], nv2[i], nv3[i], nv4[i])
                    Dx[i] <- det(cbind(nor, tetra[, 2:3], rep(1, 
                      4)))
                    Dy[i] <- -det(cbind(nor, tetra[, c(1, 3)], 
                      rep(1, 4)))
                    Dz[i] <- det(cbind(nor, tetra[, 1:2], rep(1, 
                      4)))
                    ct[i] <- det(cbind(nor, tetra[, 1:3]))
                    a[i] <- det(cbind(tetra[, 1:3], rep(1, 4)))
                  }
                  rho.sq[ind] <- 0.25 * (Dx^2 + Dy^2 + Dz^2 - 
                    4 * a * ct)/a^2
                }
                rho4 <- sqrt(rho.sq)
                rf1 <- rho4[in.tc[i1]]
                rf2 <- rho4[in.tc[i2]]
                ntri2 <- length(rf1)
                mu3 <- numeric(ntri2)
                storage.mode(mu3) <- "double"
                mu3up <- numeric(ntri2)
                storage.mode(mu3up) <- "double"
                mus3 <- .C("int3", as.integer(ntri2), as.double(rf1), 
                  as.double(rf2), mu3, mu3up, PACKAGE = "alphashape3d")
                mu3 <- c(mus3[[4]], rho4[in.tc[ih]])
                Mu3 <- c(mus3[[5]], rho4[in.tc[ih]])
                sign <- rep(c(1, -1, 1, -1), each = ntc)
                pu1 <- sign[i1] * m2340[in.tc[i1]] * m230[1:ntris] + 
                  sign[i1] * m1340[in.tc[i1]] * m130[1:ntris] + 
                  sign[i1] * m1240[in.tc[i1]] * m120[1:ntris] - 
                  2 * sign[i1] * m1230[in.tc[i1]] * m123[1:ntris]
                pu2 <- sign[i2] * m2340[in.tc[i2]] * m230[1:ntris] + 
                  sign[i2] * m1340[in.tc[i2]] * m130[1:ntris] + 
                  sign[i2] * m1240[in.tc[i2]] * m120[1:ntris] - 
                  2 * sign[i2] * m1230[in.tc[i2]] * m123[1:ntris]
                at3 <- (pu1 > 0 | pu2 > 0)
                pu3 <- sign[ih] * m2340[in.tc[ih]] * m230[(ntris + 
                  1):ntri] + sign[ih] * m1340[in.tc[ih]] * m130[(ntris + 
                  1):ntri] + sign[ih] * m1240[in.tc[ih]] * m120[(ntris + 
                  1):ntri] - 2 * sign[ih] * m1230[in.tc[ih]] * 
                  m123[(ntris + 1):ntri]
                at3 <- c(at3, pu3 > 0)
                auxmu3 <- (1 - at3) * rho3 + at3 * mu3
                storage.mode(in.trio) <- "integer"
                storage.mode(dup) <- "integer"
                storage.mode(auxmu3) <- "double"
                storage.mode(Mu3) <- "double"
                mu2 <- numeric(ned)
                storage.mode(mu2) <- "double"
                mu2up <- numeric(ned)
                storage.mode(mu2up) <- "double"
                aux1 <- 1:(3 * ntri)
                aux1 <- ceiling(aux1[ix2]/ntri)
                storage.mode(aux1) <- "integer"
                storage.mode(num.rho2) <- "double"
                storage.mode(mk0) <- "double"
                at <- numeric(ned)
                storage.mode(at) <- "integer"
                mus2 <- .C("int2", dup, as.integer(ned), as.integer(ntri), 
                  in.trio, auxmu3, Mu3, mu2, mu2up, tri.aux, 
                  aux1, ind.ed, num.rho2, mk0, at, PACKAGE = "alphashape3d")
                mu2 <- mus2[[7]]
                Mu2 <- mus2[[8]]
                at2 <- mus2[[14]]
                auxmu2 <- (1 - at2) * rho2 + at2 * mu2
                storage.mode(in.edo) <- "integer"
                storage.mode(dup1) <- "integer"
                storage.mode(auxmu2) <- "double"
                storage.mode(Mu2) <- "double"
                mu1 <- numeric(nvt)
                storage.mode(mu1) <- "double"
                mu1up <- numeric(nvt)
                storage.mode(mu1up) <- "double"
                mus1 <- .C("int1", dup1, as.integer(nvt), as.integer(ned), 
                  in.edo, auxmu2, Mu2, mu1, mu1up, PACKAGE = "alphashape3d")
                mu1 <- mus1[[7]]
                Mu1 <- mus1[[8]]
                tc.def <- cbind(tc, rho4)
                tri.def <- cbind(t1u, t2u, t3u, on.ch3, at3, 
                  rho3, mu3, Mu3)
                ed.def <- cbind(e1u, e2u, on.ch2, at2, rho2, 
                  mu2, Mu2)
                vt.def <- cbind(1:nvt, on.ch1, mu1, Mu1)
                colnames(tc.def) <- c("v1", "v2", "v3", "v4", 
                  "rhoT")
                colnames(tri.def) <- c("tr1", "tr2", "tr3", "on.ch", 
                  "attached", "rhoT", "muT", "MuT")
                colnames(ed.def) <- c("ed1", "ed2", "on.ch", 
                  "attached", "rhoT", "muT", "MuT")
                colnames(vt.def) <- c("v1", "on.ch", "muT", "MuT")
            }
            if (flag == 1) {
                flag2 <- 1
                warning("The general position assumption is not satisfied\nPerturbation of the data set was required", 
                  call. = FALSE, immediate. = TRUE)
                x <- x + rnorm(length(x), sd = sd(as.numeric(x)) * 
                  eps)
            }
        }
    }
    if (length(alphaux) > 0) {
        for (i in 1:length(alphaux)) {
            alpha <- alphaux[i]
            fclass <- numeric(length = length(tc.def[, "rhoT"]))
            fclass[alpha > tc.def[, "rhoT"]] <- 1
            tc.def <- cbind(tc.def, fclass)
            colnames(tc.def)[length(colnames(tc.def))] <- paste("fc:", 
                alpha, sep = "")
            fclass <- numeric(length = dim(tri.def)[1])
            fclass[tri.def[, "on.ch"] == 0 & tri.def[, "attached"] == 
                0 & alpha > tri.def[, "rhoT"] & alpha < tri.def[, 
                "muT"]] <- 3
            fclass[tri.def[, "on.ch"] == 0 & alpha > tri.def[, 
                "muT"] & alpha < tri.def[, "MuT"]] <- 2
            fclass[tri.def[, "on.ch"] == 0 & alpha > tri.def[, 
                "MuT"]] <- 1
            fclass[tri.def[, "on.ch"] == 1 & tri.def[, "attached"] == 
                0 & alpha > tri.def[, "rhoT"] & alpha < tri.def[, 
                "muT"]] <- 3
            fclass[tri.def[, "on.ch"] == 1 & alpha > tri.def[, 
                "muT"]] <- 2
            tri.def <- cbind(tri.def, fclass)
            colnames(tri.def)[length(colnames(tri.def))] <- paste("fc:", 
                alpha, sep = "")
            fclass <- numeric(length = dim(ed.def)[1])
            fclass[ed.def[, "on.ch"] == 0 & ed.def[, "attached"] == 
                0 & alpha > ed.def[, "rhoT"] & alpha < ed.def[, 
                "muT"]] <- 3
            fclass[ed.def[, "on.ch"] == 0 & alpha > ed.def[, 
                "muT"] & alpha < ed.def[, "MuT"]] <- 2
            fclass[ed.def[, "on.ch"] == 0 & alpha > ed.def[, 
                "MuT"]] <- 1
            fclass[ed.def[, "on.ch"] == 1 & ed.def[, "attached"] == 
                0 & alpha > ed.def[, "rhoT"] & alpha < ed.def[, 
                "muT"]] <- 3
            fclass[ed.def[, "on.ch"] == 1 & alpha > ed.def[, 
                "muT"]] <- 2
            ed.def <- cbind(ed.def, fclass)
            colnames(ed.def)[length(colnames(ed.def))] <- paste("fc:", 
                alpha, sep = "")
            fclass <- numeric(length = dim(vt.def)[1])
            fclass[alpha < vt.def[, "muT"]] <- 3
            fclass[vt.def[, "on.ch"] == 0 & alpha > vt.def[, 
                "muT"] & alpha < vt.def[, "MuT"]] <- 2
            fclass[vt.def[, "on.ch"] == 0 & alpha > vt.def[, 
                "MuT"]] <- 1
            fclass[vt.def[, "on.ch"] == 1 & alpha > vt.def[, 
                "muT"]] <- 2
            vt.def <- cbind(vt.def, fclass)
            colnames(vt.def)[length(colnames(vt.def))] <- paste("fc:", 
                alpha, sep = "")
        }
    }
    if (inh) {
        alphaux <- unique(c(alphaold, alphaux))
    }
    if (flag2 == 1) {
        ashape3d.obj <- list(tetra = tc.def, triang = tri.def, 
            edge = ed.def, vertex = vt.def, x = xorig, alpha = alphaux, 
            xpert = x)
    }
    else {
        ashape3d.obj <- list(tetra = tc.def, triang = tri.def, 
            edge = ed.def, vertex = vt.def, x = x, alpha = alphaux)
    }
    class(ashape3d.obj) <- "ashape3d"
    invisible(ashape3d.obj)
}
