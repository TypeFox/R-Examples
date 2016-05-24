SDScan <-
function (namefi = NULL, xleg = NULL, yleg = NULL, titl = NULL, 
    onecol = NULL, daty = NULL) 
{
    cat("\n*************************************************************************************\n")
    cat("\nSerial Data Scanner V1.5 \n\n")
    cat("R function for change-point detection through the Lanzante's method (Lanzante,1996)\n\n")
    cat("J. R. Lanzante, (1996). Resistant, robust and non-parametric techniques for the\n")
    cat("analysis of climate data : theory and examples, including applications to\n")
    cat("historical radiosonde station data, International Journal of Climatology,\n")
    cat("vol. 16, 1197-1226.\n")
    cat("\nOther reference : D. Amorese, (2007). Applying a change-point detection method\n")
    cat("on frequency-magnitude distributions, Bulletin of the Seismological Society of\n")
    cat("America, 97(5):1742-1749\n\n")
    cat("*************************************************************************************\n")
    NMAXITER <- 5
    SNRMIN <- 0.05
    SEUILPROB <- 0.05
    ECART <- 2
    ECARTBORD <- 2
    ECARTNV <- 0
    NLIM = 50
    FILEOUT = "SDS.res"
    no <- function(answer) answer == "n"
    yes <- function(answer) answer == "y"
    "medpairwise" <- function(n, x, y) {
        a = numeric()
        NTOTAL = 320000
        if ((n * (n - 1)/2) > NTOTAL) 
            stop("\nMedpairwise error! : too many pairs of points!\n")
        np = 0
        if (n >= 1) 
            for (i in 1:(n - 1)) {
                for (j in (i + 1):n) {
                  if (x[j] != x[i]) {
                    np = np + 1
                    if (np == NTOTAL) 
                      stop("\nToo many pairs of points...STOP!!\n")
                    a[np] = (y[j] - y[i])/(x[j] - x[i])
                  }
                }
            }
        cat("\nMEDPAIRWISE :", n, "points -> ", np, "2-points slopes\n")
        med = median(a, na.rm = T)
        b = med
        res = y - med * x
        med = median(res, na.rm = T)
        A = med
        cat("\n\tslope=", b, "intercept=", A, "\n")
        return(list(b = b, A = A))
    }
    biwmean <- function(x, c = 5, eps = 1e-04) {
        m <- median(x)
        s <- median(abs(x - m))
        u <- (x - m)/(c * s + eps)
        w <- rep(0, length(x))
        i <- abs(u) <= 1
        w[i] <- ((1 - u^2)^2)[i]
        bm <- sum(w * x)/sum(w)
        return(bm)
    }
    biwvar <- function(x, c = 5, eps = 1e-04) {
        m <- median(x)
        s <- median(abs(x - m))
        u <- (x - m)/(c * s + eps)
        w <- rep(0, length(x))
        i <- abs(u) <= 1
        w[i] <- (1 - u^2)[i]
        w5 <- 1 - 5 * u^2
        bwm <- sqrt(length(x)) * sqrt(sum(w^4 * (x - m)^2))/abs(sum(w * 
            w5))
        return(bwm)
    }
    "quickregli" <- function(n, x, y) {
        sumx = sum(x)
        sumx2 = sum(x * x)
        sumy = sum(y)
        sumy2 = sum(y * y)
        sumxy = sum(x * y)
        delta = n * sumx2 - sumx * sumx
        b = (n * sumxy - sumx * sumy)/delta
        A = (sumy * sumx2 - sumx * sumxy)/delta
        return(list(b = b, A = A))
    }
    "ranktest" <- function(niter, N, yvar) {
        SA = numeric()
        corsa = (N + 1) * seq(1:N)
        SA = abs(2 * cumsum(rank(yvar, ties.method = "average")) - 
            corsa)
        n1 <- max(which(SA == unique(SA)[order(unique(SA))[length(unique(SA)) - 
            niter + 1]]))
        W = sum(rank(yvar, ties.method = "average")[1:n1])
        cat("\nn1=", n1, "W=", W, "\n")
        n2 = N - n1
        Wcrit = n1 * (N + 1)/2
        cat("\n niter=", niter)
        cat("n=", N, "n2=", n2, "Critical W=", Wcrit, "\n")
        sw = sqrt(n1 * n2 * (N + 1)/12)
        srang2 = sum((rank(yvar, ties.method = "average"))^2)
        srang_star = sum(rank(yvar, ties.method = "first")[n1:N])
        sw2 = sqrt((n1 * n2 * srang2)/(N * (N - 1)) - (n1 * n2 * 
            (N + 1) * (N + 1))/(4 * (N - 1)))
        if (n1 <= n2) 
            wx = W
        if (n2 < n1) 
            wx = srang_star
        xn1 <- yvar[1:n1]
        xn2 <- yvar[-(1:n1)]
        wt = wilcox.test(xn1, xn2, exact = F, correct = T)
        p = wt$p.value
        cat("sw=", sw, "sw_t=", sw2, "W=", W, "Wx=", wx, "\n")
        return(list(p = p, tau = n1))
    }
    "snr" <- function(ndat, ndis, t, x, w, tabn, choix, NLIM) {
        ntabn = numeric()
        w1 = numeric()
        w2 = numeric()
        w3 = numeric()
        x1 = numeric()
        x2 = numeric()
        x3 = numeric()
        rdn = 1
        ntabn = tabn
        ntabn = c(ntabn, t)
        ii = max(rank(ntabn, ties.method = "first")[ntabn == 
            t])
        nl = ntabn[rank(ntabn, ties.method = "first") == (ii - 
            1)]
        if (nl == t && nl != ntabn[1]) 
            nl = ntabn[rank(ntabn, ties.method = "average") == 
                (ii - 2)]
        nr = ntabn[rank(ntabn, ties.method = "average") == (ii + 
            1)]
        if (nr == t && nr != ntabn[2]) 
            nr = ntabn[rank(ntabn, ties.method = "average") == 
                (ii + 2)]
        if (nl == 1) 
            N1 = t - nl + 1
        if (nl != 1) 
            N1 = t - nl
        N2 = nr - t
        if (N1 > NLIM) {
            N1 = NLIM
            nl = t - N1
        }
        if (N2 > NLIM) {
            N2 = NLIM
            nr = t + N2
        }
        n = N1 + N2
        if (nl == 1) {
            x1[1:N1] = x[1:N1]
            w1[1:N1] = w[1:N1]
        }
        if (nl != 1) {
            x1[1:N1] = x[(nl + 1):(nl + N1)]
            w1[1:N1] = w[(nl + 1):(nl + N1)]
        }
        x2[1:N2] = x[(t + 1):(t + N2)]
        w2[1:N2] = w[(t + 1):(t + N2)]
        x3 = c(x1, x2)
        w3 = c(w1, w2)
        if (n <= 800) {
            mpw = medpairwise(n, w3, x3)
            b = mpw$b
            A = mpw$A
        }
        if (n > 800) {
            qrl = quickregli(n, w3, x3)
            b = qrl$b
            A = qrl$A
        }
        if (choix == 1) 
            cat("\nb=", b, "A=", A, "\n")
        xl = biwmean(x1, 7.5, 1e-04)
        xr = biwmean(x2, 7.5, 1e-04)
        xleft = xl
        xright = xr
        omean = (N1 * xl + N2 * xr)/n
        var = (N1 * (xl - omean) * (xl - omean) + N2 * (xr - 
            omean) * (xr - omean))/(n - 1)
        if (choix != 1) 
            x1[1:N1] = x1[1:N1] - xl
        if (choix == 1) 
            x1[1:N1] = x1[1:N1] - b * (w1 - w[t])
        if (choix != 1) 
            x2[1:N2] = x2[1:N2] - xr
        if (choix == 1) 
            x2[1:N2] = x2[1:N2] - b * (w2 - w[t])
        x3 = c(x1, x2)
        noisesd = biwvar(x3, 7.5, 1e-04)
        noisev = noisesd * noisesd
        if (noisev == 0) 
            noisev = 1e-06
        if (choix == 0 || choix == 1) 
            rdn = noisev
        if (choix == 2) 
            rdn = var/noisev
        cat("\nsnr -> xl=", xl, "n=", N1, "nl=", nl, "xr=", xr, 
            "n=", N2, "nr=", nr, "x=", omean, "sd=", var, "sn=", 
            noisev, "snr=", rdn, "\n")
        ratio = rdn
        return(list(xleft = xleft, xright = xright, ratio = ratio))
    }
    "dtrend" <- function(ndat, ndis, t, x, w, tabn, NLIM) {
        ntabn = numeric()
        X = numeric()
        x3 = numeric()
        w2 = numeric()
        ntabn[1:ndis] = tabn[1:ndis]
        ntabn = c(ntabn, t)
        ii = max(rank(ntabn, ties.method = "average")[ntabn == 
            t])
        nl = ntabn[rank(ntabn, ties.method = "average") == (ii - 
            1)]
        if (nl == t && nl != ntabn[1]) 
            nl = ntabn[rank(ntabn, ties.method = "average") == 
                (ii - 2)]
        nr = ntabn[rank(ntabn, ties.method = "average") == (ii + 
            1)]
        if (nr == t && nr != ntabn[2]) 
            nr = ntabn[rank(ntabn, ties.method = "average") == 
                (ii + 2)]
        if (nl == 1) 
            N1 = t - nl + 1
        if (nl != 1) 
            N1 = t - nl
        N2 = nr - t
        if (N1 > NLIM) {
            N1 = NLIM
            nl = t - N1
        }
        if (N2 > NLIM) {
            N2 = NLIM
            nr = t + N2
        }
        n = N1 + N2
        if (nl == 1) {
            x3[1:n] = x[(nl):(nl + n - 1)]
            w2[1:n] = w[(nl):(nl + n - 1)]
        }
        if (nl != 1) {
            x3[1:n] = x[(nl + 1):(nl + n)]
            w2[1:n] = w[(nl + 1):(nl + n)]
        }
        if (n <= 800) {
            mpw = medpairwise(n, w2, x3)
            b = mpw$b
            A = mpw$A
        }
        if (n > 800) {
            qrl = quickregli(n, w2, x3)
            b = qrl$b
            A = qrl$A
        }
        cat("\nReduction :", b, A, "\n")
        X[1:nl] = x[1:nl]
        X[(nl + 1):nr] = x[(nl + 1):nr] + b * (t - w[i])
        X[(nr + 1):ndat] = x[(nr + 1):ndat]
        x3[1:ndat] = X[1:ndat]
        return(list(dt = x3))
    }
    "adjust" <- function(dat, m, tabn) {
        ntabn = numeric()
        idebut = numeric()
        ifin = numeric()
        mediane = numeric()
        jdeb = numeric()
        jfin = numeric()
        ntabn[1:m] = tabn[1:m]
        tabsort = sort(ntabn)
        for (i in 1:(m - 1)) {
            dat2 = numeric()
            if (i == 1) 
                dat2 = dat[tabsort[i]:tabsort[i + 1]]
            if (i != 1) 
                dat2 = dat[(tabsort[i] + 1):tabsort[i + 1]]
            med = median(dat2, na.rm = T)
            if (i == 1) {
                jdeb[i] = tabsort[i]
                jfin[i] = tabsort[i + 1]
                for (j in jdeb[i]:jfin[i]) {
                  dat[j] = dat[j] - med
                }
            }
            if (i != 1) {
                jdeb[i] = 1 + tabsort[i]
                jfin[i] = tabsort[i + 1]
                for (j in jdeb[i]:jfin[i]) {
                  dat[j] = dat[j] - med
                }
            }
            idebut[i] = jdeb[i]
            ifin[i] = jfin[i]
            mediane[i] = med
        }
        return(list(ndat = dat, debut = idebut, fin = ifin, mediane = mediane))
    }
    readlineB <- function(myString0) {
        print(myString0)
        a <- scan(what = character(0), nmax = 1)
        return(a)
    }
    if (is.null(namefi)) {
        namefi = readlineB("\nData file name ? : ")
    }
    b = read.table(namefi)
    if (is.null(xleg)) {
        xleg = "X"
        xleg = readlineB("\nX-axis title ? : ")
    }
    if (is.null(yleg)) {
        yleg = "Y"
        yleg = readlineB("\nY-axis title ? : ")
    }
    if (is.null(titl)) {
        titl = "SERIES"
        titl = readlineB("\nMain title ? : ")
    }
    if (is.null(onecol)) {
        onecol = "n"
        onecol = readlineB("\nIs this a single column file ? (y/n): ")
    }
    if (no(onecol)) 
        xvar <- b$V1
    if (yes(onecol)) 
        xvar <- 1:length(b$V1)
    if (is.null(daty)) {
        daty = "n"
        daty = readlineB("\nShould gradient values be processed ? (y/n): ")
    }
    if (yes(onecol)) 
        b$V2 <- b$V1
    if (no(daty)) 
        yvar <- b$V2
    if (yes(daty)) {
        yvar <- diff(b$V2)/diff(b$V1)
        xvar <- xvar[-length(xvar)]
        xvar <- xvar[yvar != Inf]
        yvar <- yvar[yvar != Inf]
    }
    pva = numeric()
    tau = numeric()
    N <- length(yvar)
    chpt = numeric()
    testchpt = numeric()
    tnvd = numeric()
    tnvt = numeric()
    tprobrob = numeric()
    maxtmp = numeric()
    cat("\nChange-point detection is being performed\n")
    for (i in 1:N) {
        chpt[i] = 0
        testchpt[i] = 0
    }
    yvar0 <- yvar
    chpt = rep(0, length(yvar))
    testchpt = rep(0, length(yvar))
    k = 1
    chpt[k] = 1
    k = 2
    chpt[k] = N
    ca <- 0
    cort <- 0
    niter <- 0
    while (niter < NMAXITER) {
        niter <- niter + 1
        cat("\nITERATION #", niter, "\n\n")
        rkt = ranktest(niter, N, yvar)
        p = rkt$p
        tau = rkt$tau
        cat("p-value =", p, "ntau_i =", tau, "\n")
        if (p < SEUILPROB) {
            cat("\tTest Passed\n")
            nxrob = tau
            nyrob = N - tau
            xn1 <- yvar[1:nxrob]
            xn2 <- yvar[-(1:nxrob)]
            z = c(xn1, xn2)
            ix = rep(0, length(xn1))
            for (i in 1:length(xn1)) {
                for (j in (length(xn1) + 1):length(z)) {
                  if (rank(z, ties.method = "average")[j] < rank(z, 
                    ties.method = "average")[i]) 
                    ix[i] = ix[i] + 1
                }
            }
            iy = rep(0, length(xn2))
            for (i in 1:length(xn2)) {
                for (j in 1:length(xn1)) {
                  if (rank(z, ties.method = "average")[j] < rank(z, 
                    ties.method = "average")[i + length(xn1)]) 
                    iy[i] = iy[i] + 1
                }
            }
            sumnx = sum(ix)
            sumnx2 = sum(ix * ix)
            sumny = sum(iy)
            sumny2 = sum(iy * iy)
            nxbar = sumnx/length(xn1)
            nybar = sumny/length(xn2)
            sdnx = sumnx2 - length(xn1) * nxbar * nxbar
            sdny = sumny2 - length(xn2) * nybar * nybar
            zstat = 0.5 * (length(xn1) * nxbar - length(xn2) * 
                nybar)/sqrt(nxbar * nybar + sdnx + sdny)
            zval = zstat
            if (zstat > 0) 
                zval = -zval
            p = pnorm(zval)
            p = p * 2
            cat("\nnxbar=", nxbar, "nybar=", nybar, "sdnx=", 
                sdnx, "sdny=", sdny, "zstat=", zstat, "p=", p, 
                "\n")
            if (length(xn1) < 12 && length(xn2) < 12) 
                p = -p
            probrob = p
            pasymp = 1
            if (probrob < 0) {
                cat("\n\tWARNING : exact p-value against asymptotic!\n")
                probrob = -probrob
                pasymp = 0
            }
            if (probrob < SEUILPROB) 
                cat("\tRobust Rank Order Test Passed\n")
            snr1 = snr(N, k, tau, yvar0, xvar, chpt, 0, NLIM)
            NVD = snr1$ratio
            snr2 = snr(N, k, tau, yvar0, xvar, chpt, 1, NLIM)
            xl = snr2$xleft
            xr = snr2$xright
            NVT = snr2$ratio
            cat("\nDisc. Noise Var.=", NVD, "Trend Noise Var=", 
                NVT, "\n")
            prox = 1
            cat("\n", tau, "<->")
            for (i in 1:k) {
                cat(" ", chpt[i])
                if (abs(tau - chpt[i]) <= ECART) 
                  prox = 0
            }
            if (prox == 1) 
                cat("\nGap between change-points")
            if ((tau - 1) <= ECARTBORD) 
                prox = 0
            if ((N - tau) <= ECARTBORD) 
                prox = 0
            if ((NVD - NVT) > ECARTNV && prox == 1) {
                cat("\tTrend reduction\n")
                cort = cort + 1
                if (cort >= NMAXITER) 
                  break
                suptend = dtrend(N, k, tau, yvar0, xvar, chpt, 
                  NLIM)
                yvar = suptend$dt
                niter = 0
            }
            if ((NVD - NVT) <= ECARTNV && prox == 1) {
                k = k + 1
                if (k > 3) {
                  for (i in 1:(k - 2)) yvar[debut[i]:fin[i]] = yvar[debut[i]:fin[i]] + 
                    medi[i]
                }
                chpt[k] = tau
                testchpt[tau] = testchpt[tau] + 1
                tnvd[k] = NVD
                tnvt[k] = NVT
                if (pasymp == 1) 
                  tprobrob[k] = probrob
                if (pasymp == 0) 
                  tprobrob[k] = -probrob
                cat("\nCHANGE-POINT DETECTED\n")
                ajt = adjust(yvar, k, chpt)
                yvar = ajt$ndat
                debut = ajt$debut
                fin = ajt$fin
                medi = ajt$mediane
                niter = 0
                ca = ca + 1
                cat("\nADJUSTMENT #", ca)
            }
        }
    }
    cat("\nNUMBER OF ITERATIONS :", niter)
    if (niter != NMAXITER) 
        if ((NVD - NVT) <= ECARTNV && prox == 1) 
            for (i in 1:(k - 1)) for (j in debut[i]:fin[i]) yvar[j] = yvar[j] + 
                medi[i]
    xseg = numeric()
    yseg = numeric()
    yseg2 = numeric()
    pseg = numeric()
    s2nr = numeric()
    ligne = character()
    for (i in 1:k) {
        snrstar = ""
        nvstar0 = ""
        pstar = ""
        nvstar = ""
        tau = chpt[i]
        if (tau != N && tau != 1) {
            snr3 = snr(N, k, tau, yvar0, xvar, chpt, 0, N)
            XL = snr3$xleft
            XR = snr3$xright
            NVD = snr3$ratio
            snr4 = snr(N, k, tau, yvar0, xvar, chpt, 1, N)
            NVT = snr4$ratio
            snr5 = snr(N, k, tau, yvar0, xvar, chpt, 2, N)
            SNRD = snr5$ratio
            if (SNRD < SNRMIN) 
                snrstar = paste(format(SNRD, digits = 6), "*", 
                  sep = "")
            if (SNRD >= SNRMIN) 
                snrstar = paste(snrstar, format(SNRD, digits = 6), 
                  sep = "")
            if ((tnvd[i] - tnvt[i]) > ECARTNV) 
                nvstar0 = paste(format(tnvt[i], digits = 6), 
                  "*", sep = "")
            if ((tnvd[i] - tnvt[i]) <= ECARTNV) 
                nvstar0 = paste(nvstar0, format(tnvt[i], digits = 6), 
                  sep = "")
            if (tprobrob[i] < 0 || tprobrob[i] >= SEUILPROB) 
                pstar = paste(format(chpt[i], digits = 6), "*", 
                  sep = "")
            if (tprobrob[i] >= 0 && tprobrob[i] < SEUILPROB) 
                pstar = paste(pstar, format(chpt[i], digits = 6), 
                  sep = "")
            if ((NVD - NVT) > ECARTNV) 
                nvstar = paste(format(NVT, digits = 6), "*", 
                  sep = "")
            if ((NVD - NVT) <= ECARTNV) 
                nvstar = paste(nvstar, format(NVT, digits = 6), 
                  sep = "")
            ligne[i] = paste(format(xvar[tau], digits = 6), " & ", 
                format(yvar0[tau], digits = 6), " & ", pstar, 
                " & ", format(tnvd[i], digits = 6), " & ", nvstar0, 
                " & ", format(NVD, digits = 6), " & ", nvstar, 
                " & ", snrstar, " & ", format(XL, digits = 6), 
                " & ", format(XR, digits = 6), sep = "")
            xseg = c(xseg, xvar[tau])
            yseg = c(yseg, XL)
            yseg2 = c(yseg2, XR)
            pseg = c(pseg, tprobrob[i])
            s2nr = c(s2nr, SNRD)
        }
    }
    nn = order(xseg)
    nn = c(nn, 0)
    xseg1 = c(min(xvar), sort(xseg))
    xseg2 = c(sort(xseg), max(xvar))
    yseg = yseg[order(xseg)]
    yseg2 = yseg2[order(xseg)]
    yseg = c(yseg, yseg2[length(yseg2)])
    pseg = pseg[order(xseg)]
    pseg = c(pseg, 0)
    s2nr = s2nr[order(xseg)]
    s2nr = c(s2nr, 0)
    dfr = data.frame(cbind(nn, xseg1, yseg, xseg2, pseg, s2nr))
    ligne[2] = "X_value & Y_value & Chpt & NoiVarDi & NoiVarTr & NVDpost & NVTpost & DisSNR & LeftBWM & RightBWM"
    ligne = ligne[!is.na(ligne)]
    write.table(ligne, FILEOUT, quote = F, col.names = F, row.names = F)
    cat("\nNumerical results in SDS.res\n\n")
    "niceplot" <- function(df, lab1, lab2, mai, df2, locleg) {
        if (substr(mai, 1, 5) == "paste") 
            mai = parse(text = mai)
        if (substr(lab1, 1, 5) == "paste") 
            lab1 = parse(text = lab1)
        if (substr(lab2, 1, 5) == "paste") 
            lab2 = parse(text = lab2)
        par(mfrow = c(1, 1))
        par(bg = "lightgray")
        par(mar = c(5, 5, 4, 2) + 0.1, xpd = TRUE)
        plot(df, type = "n", axes = FALSE, ann = FALSE)
        usr = par("usr")
        rect(usr[1], usr[3], usr[2], usr[4], col = "cornsilk", 
            border = "black")
        lines(df, col = "blue")
        points(df, pch = 21, bg = "lightcyan", cex = 1.25)
        if (length(df2[, 2]) > 1) {
            segments(df2[, 2], df2[, 3], df2[, 4], lwd = 4)
            ntext = length(df2[, 2])
            text(df2[-ntext, 4], df2[-ntext, 3], as.character(df2[-ntext, 
                1]), col = "red", cex = 1.7, pos = 4, offset = 0.35)
        }
        axis(2, col.axis = "blue", las = 1)
        axis(1, col.axis = "blue")
        box()
        if (length(df2[, 2]) > 1) {
            nline = rep("", ntext)
            if (locleg[1] == 1) {
                cat("\nPLEASE, locate with the mouse the topright corner of the legend in the plot window\n\n")
                loc = locator(1)
                if (!is.null(loc)) 
                  outloc = c(loc$x, loc$y)
                if (is.null(loc)) {
                  loc = "topleft"
                  outloc = loc
                }
            }
            if (locleg[1] != 1) 
                loc = locleg
            temp <- legend(loc[1], loc[2], inset = c(0, 0), legend = nline, 
                xjust = 1, yjust = 1, title = "    Statistics for change-points    ", 
                cex = 0.8)
            str = sprintf("%.2e", df2[-ntext, 5])
            nch = nchar(str)
            textline = sprintf("%2d %9.2f %7.2f %s", df2[-ntext, 
                1], df2[-ntext, 4], df2[-ntext, 6], format(df2[-ntext, 
                5], digits = 11 - nch, width = 9))
            hdr = sprintf("%2s %9s %7s %9s", "N", "XChpt", "SNR", 
                "P-value")
            textline = c(hdr, textline)
            text(temp$rect$left + temp$rect$w, temp$text$y, textline, 
                pos = 2, cex = 0.8)
        }
        title(main = mai, font.main = 4, col.main = "red", cex.main = 1.7)
        title(xlab = lab1, col.lab = "red", cex.lab = 1.4)
        title(ylab = lab2, col.lab = "red", cex.lab = 1.4)
        if (length(df2[, 2]) <= 1) 
            outloc = c("No detection")
        if (locleg[1] == 1) 
            return(list(inloc = outloc))
    }
    pts = data.frame(xvar, yvar0)
    nplt = niceplot(pts, xleg, yleg, titl, dfr, 1)
    inloc = nplt$inloc
    cat(inloc)
    resnamepdf = "SDS.pdf"
    pdf(resnamepdf, version = "1.4")
    niceplot(pts, xleg, yleg, titl, dfr, inloc)
    dev.off()
    resnamepng = "SDS.png"
    png(resnamepng, width = 1200, height = 1200, res = 120)
    niceplot(pts, xleg, yleg, titl, dfr, inloc)
    dev.off()
    cat("\nGraphics in SDS.png\n\n")
    cat("\nGraphics in SDS.pdf\n\n")
}
