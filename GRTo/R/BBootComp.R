BBootComp <-
function (finame1, finame2 = NULL, colid1 = 1, colid2 = 1, hd1 = FALSE, 
    hd2 = FALSE, nrep = 5000, alter = c("two.sided", "less", 
        "greater"), tm1 = NULL, tm2 = NULL, findtm1 = TRUE, findtm2 = TRUE, 
    plot = FALSE, title = "BootComp", oplt = "p") 
{
    alter <- match.arg(alter)
    rd = 1
    bw = 0.1
    "m0b" <- function(a, delta = 0.1, tm = NULL, findtm = TRUE) {
        tau <- numeric()
        pva <- numeric()
        minmag <- min(a, na.rm = TRUE)
        maxmag <- max(a, na.rm = TRUE)
        g_r <- hist(a, plot = FALSE, breaks = round(seq((minmag - 
            delta/2), (maxmag + delta/2), delta), 3))
        n <- length(g_r$density)
        xc <- round(seq(minmag, max(a, na.rm = TRUE), delta)[1:(n - 
            1)], 3)
        x <- round(seq(minmag, max(a, na.rm = TRUE), delta), 
            3)
        log_n <- log10((1/delta) * g_r$counts * delta)
        x <- x[is.finite(log_n)]
        log_n <- log_n[is.finite(log_n)]
        sl <- diff(log_n)/diff(x)
        xsl <- x[2:length(x)]
        niter <- 3
        N <- length(sl)
        j <- 0
        k <- 0
        SA <- vector(length = N)
        while (j < niter) {
            for (i in seq(1, N, 1)) SA[i] <- abs(2 * sum(rank(sl)[1:i]) - 
                i * (N + 1))
            n1 <- which(SA == SA[order(SA)[length(order(SA))]])
            xn1 <- sl[1:n1[1]]
            xn2 <- sl[-(1:n1[1])]
            if ((n1[1] > 2) & (n1[1] <= (N - 2)) & (wilcox.test(xn1, 
                xn2, exact = FALSE, correct = TRUE)[3] < 0.05)) {
                k <- k + 1
                pva[k] <- wilcox.test(xn1, xn2, exact = FALSE, 
                  correct = TRUE)[3]
                tau[k] <- n1[1]
                if (k > 1) {
                  medsl1 <- median(sl[1:n0])
                  medsl2 <- median(sl[-(1:n0)])
                  for (i in seq(1, n0, 1)) sl[i] <- sl[i] + medsl1
                  for (i in seq(n0 + 1, length(sl), 1)) sl[i] <- sl[i] + 
                    medsl2
                }
                medsl1 <- median(sl[1:n1[1]])
                medsl2 <- median(sl[-(1:n1[1])])
                for (i in seq(1, n1[1], 1)) sl[i] <- sl[i] - 
                  medsl1
                for (i in seq(n1[1] + 1, length(sl), 1)) sl[i] <- sl[i] - 
                  medsl2
                n0 <- n1[1]
            }
            j <- j + 1
        }
        v_pva = as.vector(pva, mode = "numeric")
        ip = order(v_pva)
        m0 = c()
        bv = c()
        m0 = c(signif(xsl[tau[ip[1]]]), signif(xsl[tau[ip[2]]]))
        if (!is.null(tm)) {
            m0 = c(tm, tm)
        }
        if (is.null(tm) & !findtm) {
            m0 = c(min(a, na.rm = TRUE), min(a, na.rm = TRUE))
        }
        if (!is.na(m0[1])) 
            bau = 1/(log(10) * (mean(a[which(a >= m0[1])], na.rm = TRUE) - 
                (m0[1] - delta/2)))
        if (is.null(tm) & is.na(m0[1]) & findtm) 
            bau = NA
        res = list()
        res$m0 = m0[1]
        res$bv = bau
        res$aval = log10(sum(a >= m0[1])) + bau * m0[1]
        res
    }
    mrand.test <- function(x, y, bw = 0.1, nrep = 5000, alternative = c("two.sided", 
        "less", "greater"), tm1 = NULL, tm2 = NULL, findtm1 = TRUE, 
        findtm2 = TRUE, save = FALSE) {
        alternative <- match.arg(alternative)
        pooled = c(x, y)
        n = length(pooled)
        nv1 = length(x)
        nv2 = length(y)
        cat("\nNumber of values (file #1, file #2) :", nv1, nv2, 
            "\n")
        cat("\n...you just have to wait...\n")
        ar = array(0, c(nrep, n))
        for (i in 1:nrep) {
            ar[i, ] = sample(pooled, n, replace = TRUE)
        }
        mx = apply(ar[, 1:length(x)], 1, m0b, delta = bw, tm = tm1, 
            findtm = findtm1)
        my = apply(ar[, (length(x) + 1):n], 1, m0b, delta = bw, 
            tm = tm2, findtm = findtm2)
        rep1 = c()
        rep12 = c()
        rep2 = c()
        rep22 = c()
        for (i in 1:nrep) {
            rep1[i] = mx[[i]]$m0[1]
            rep12[i] = mx[[i]]$bv
            rep2[i] = my[[i]]$m0[1]
            rep22[i] = my[[i]]$bv
        }
        cat("\nNumber of valid replicates (null hypothesis=pooled magnitude values):\n")
        cat("Synthetic set 1 :", nrep - sum(is.na(rep12)))
        cat("\nSynthetic set 2 :", nrep - sum(is.na(rep22)))
        thetatmp = rep12 - rep22
        theta = thetatmp[!is.na(thetatmp)]
        nrep = length(theta)
        cat("\nTotal :", nrep, "\n")
        cat("\nQuantiles of the results for the bootstrap sample \n(null hypothesis=pooled magnitude values):\n")
        cat("Threshold magnitude values:\n")
        cat(signif(quantile(rep1, probs = c(0.025, 0.5, 0.975), 
            na.rm = TRUE), 2), "\n")
        cat(signif(quantile(rep2, probs = c(0.025, 0.5, 0.975), 
            na.rm = TRUE), 2), "\n\n")
        cat("b-values:\n")
        cat(signif(quantile(rep12, probs = c(0.025, 0.5, 0.975), 
            na.rm = TRUE), 2), "\n")
        cat(signif(quantile(rep22, probs = c(0.025, 0.5, 0.975), 
            na.rm = TRUE), 2), "\n")
        cat("\n...you just have to wait...\n")
        ar1 = array(0, c(nrep, nv1))
        for (i in 1:nrep) {
            ar1[i, ] = sample(x, nv1, replace = TRUE)
        }
        ar2 = array(0, c(nrep, nv2))
        for (i in 1:nrep) {
            ar2[i, ] = sample(y, nv2, replace = TRUE)
        }
        mv1 = apply(ar1, 1, m0b, delta = bw, tm = tm1, findtm = findtm1)
        mv2 = apply(ar2, 1, m0b, delta = bw, tm = tm2, findtm = findtm2)
        mag1 = c()
        bv1 = c()
        av1 = c()
        mag2 = c()
        bv2 = c()
        av2 = c()
        for (i in 1:nrep) {
            mag1[i] = mv1[[i]]$m0[1]
            bv1[i] = mv1[[i]]$bv
            av1[i] = mv1[[i]]$aval
            mag2[i] = mv2[[i]]$m0[1]
            bv2[i] = mv2[[i]]$bv
            av2[i] = mv2[[i]]$aval
        }
        cat("\nNumber of valid replicates (observed values):\n")
        cat("Observed set 1 :", nrep - sum(is.na(bv1)))
        cat("\nObserved set 2 :", nrep - sum(is.na(bv2)), "\n")
        if (!is.null(tm1)) 
            cat("\nFixed threshold magnitude value for set 1 :", 
                tm1)
        if (is.null(tm1)) {
            cat("\nQuantiles for threshold magnitude values, observed set 1 :\n")
            cat(signif(quantile(mag1, probs = c(0.025, 0.5, 0.975), 
                na.rm = TRUE), 2))
            tm1 = quantile(mag1, probs = 0.5, na.rm = TRUE)
        }
        cat("\nNumber of events above m_0 for set 1 :", (sum(x >= 
            tm1)), "\n")
        if (!is.null(tm2)) 
            cat("\nFixed threshold magnitude value for set 2 :", 
                tm2)
        if (is.null(tm2)) {
            cat("\nQuantiles for threshold magnitude values, observed set 2 :\n")
            cat(signif(quantile(mag2, probs = c(0.025, 0.5, 0.975), 
                na.rm = TRUE), 2))
            tm2 = quantile(mag2, probs = 0.5, na.rm = TRUE)
        }
        cat("\nNumber of events above m_0 for set 2 :", (sum(y >= 
            tm2)))
        cat("\n\nQuantiles for b-values, observed set 1 :\n")
        cat(signif(quantile(bv1, probs = c(0.025, 0.5, 0.975), 
            na.rm = TRUE), 2))
        cat("\n\nQuantiles for b-values, observed set 2 :\n")
        cat(signif(quantile(bv2, probs = c(0.025, 0.5, 0.975), 
            na.rm = TRUE), 2))
        b1 = median(bv1, na.rm = TRUE)
        b2 = median(bv2, na.rm = TRUE)
        sd1 = sd(bv1, na.rm = TRUE)
        sd2 = sd(bv2, na.rm = TRUE)
        a1 = median(av1, na.rm = TRUE)
        a2 = median(av2, na.rm = TRUE)
        bvar = sd1^2 + sd2^2
        cat("\n\nStandard error of bootstrap sample for the b-value \nobserved set 1 :", 
            sd1)
        cat("\n\nStandard error of bootstrap sample for the b-value \nobserved set 2 :", 
            sd2)
        cat("\n\nVariance of the bootstrap sample difference :", 
            bvar, "\n")
        zval = (b1 - b2)/sqrt(bvar)
        cat("***\n")
        tobs = b1 - b2
        cat("The true number of replicates (probability computations) is", 
            nrep, "\n\n")
        p1 = sum(theta >= tobs)/nrep
        p2 = 1 - p1
        p3 = sum(abs(theta) >= abs(tobs), na.rm = TRUE)/nrep
        p = switch(alternative, greater = p1, less = p2, two.sided = p3)
        mt = mean(theta, na.rm = TRUE)
        ci = quantile(theta, c(0.025, 0.975), na.rm = TRUE)
        if (save) {
            write.table(signif(theta, 4), file = "replic.res", 
                row.names = FALSE, col.names = FALSE)
            cat("\nfile replic.rec created (bootstrap differences)\n")
        }
        names(tobs) <- "b-value difference"
        val = list(statistic = tobs, parameter = NULL, p.value = as.numeric(p), 
            alternative = alternative, method = "bootstrap hypothesis testing", 
            data.name = paste(finame1, "and", finame2))
        class(val) <- "htest"
        return(list(val = val, b1 = b1, b2 = b2, sd1 = sd1, sd2 = sd2, 
            m01 = tm1, m02 = tm2, a1 = a1, a2 = a2))
    }
    brand <- function(x, bw = 0.1, nrep = 200, tm = NULL, findtm = TRUE) {
        nv1 = length(x)
        cat("\nNumber of events :", nv1, "\n")
        cat("\n...you just have to wait...\n")
        ar1 = array(0, c(nrep, nv1))
        for (i in 1:nrep) ar1[i, ] = sample(x, nv1, replace = TRUE)
        mv1 = apply(ar1, 1, m0b, delta = bw, tm = tm, findtm = findtm)
        mag1 = c()
        bv1 = c()
        av1 = c()
        for (i in 1:nrep) {
            mag1[i] = mv1[[i]]$m0[1]
            bv1[i] = mv1[[i]]$bv
            av1[i] = mv1[[i]]$aval
        }
        cat("\nNumber of valid replicates :", nrep - sum(is.na(bv1)), 
            "/", nrep, "\n", sep = "")
        cat("These replicates concern computations of m_0")
        cat("\nQuantiles for threshold magnitude values :\n")
        cat(signif(quantile(mag1, probs = c(0.025, 0.5, 0.975), 
            na.rm = TRUE), 2))
        tm = quantile(mag1, probs = 0.5, na.rm = TRUE)
        cat("\nNumber of events above m_0 :", (sum(x >= tm)))
        cat("\nQuantiles for b-values :\n")
        cat(signif(quantile(bv1, probs = c(0.025, 0.5, 0.975), 
            na.rm = TRUE), 2))
        b = median(bv1, na.rm = TRUE)
        sd = sd(bv1, na.rm = TRUE)
        a = median(av1, na.rm = TRUE)
        cat("\nMedian and standard error of bootstrap sample for the b-value :\n")
        return(list(b = b, sd = sd, m0 = tm, a = a))
    }
    m1 = read.table(file = finame1, header = hd1)
    if (!is.null(finame2)) 
        m2 = read.table(file = finame2, header = hd2)
    mtmp = m1[, colid1]
    m1 = mtmp[!is.na(mtmp)]
    if (rd != 0) 
        m1 = round(m1, rd)
    if (any(m1 <= -3) | any(m1 >= 9.5)) 
        stop("\nAre you sure that input values are magnitude values ? Wrong input column ?\n\n")
    if (!is.null(finame2)) {
        mtmp = m2[, colid2]
        m2 = mtmp[!is.na(mtmp)]
        if (rd != 0) 
            m2 = round(m2, rd)
        if (any(m2 <= -3) | any(m2 >= 9.5)) 
            stop("\nAre you sure that input values are magnitude values ? Wrong input column ?\n\n")
    }
    if (!is.null(finame2)) {
        randomisation = mrand.test(m1, m2, bw, nrep, alternative = alter, 
            tm1, tm2, findtm1, findtm2)
        b1 = randomisation$b1
        sd1 = randomisation$sd1
        b2 = randomisation$b2
        sd2 = randomisation$sd2
        m01 = randomisation$m01
        m02 = randomisation$m02
    }
    if (is.null(finame2)) {
        randomisation = brand(m1, bw, nrep, tm1, findtm1)
        b1 = randomisation$b
        sd1 = randomisation$sd
        m01 = randomisation$m0
    }
    print(randomisation)
    if (plot) {
        minmag <- min(m1, na.rm = TRUE)
        maxmag <- max(m1, na.rm = TRUE)
        g_r <- hist(m1, plot = FALSE, breaks = round(seq((minmag - 
            bw/2), (maxmag + bw/2), bw), 3))
        n <- length(g_r$density)
        xc1 <- round(seq(minmag, max(m1, na.rm = TRUE), bw)[1:(n - 
            1)], 3)
        y1 = numeric()
        for (i in 1:length(xc1)) y1[i] = sum(m1 >= xc1[i])
        if (!is.null(finame2)) {
            minmag <- min(m2, na.rm = TRUE)
            maxmag <- max(m2, na.rm = TRUE)
            g_r <- hist(m2, plot = FALSE, breaks = round(seq((minmag - 
                bw/2), (maxmag + bw/2), bw), 3))
            n <- length(g_r$density)
            xc2 <- round(seq(minmag, max(m2, na.rm = TRUE), bw)[1:(n - 
                1)], 3)
            y2 = numeric()
            for (i in 1:length(xc2)) y2[i] = sum(m2 >= xc2[i])
        }
        figname = paste(title, "_fig.png", sep = "")
        png(filename = figname, width = 2400, height = 2715, 
            res = 360)
        ind1 = which(round(xc1, 6) == m01)
        A1 = log10(y1[ind1]) + b1 * m01
        if (!is.null(finame2)) {
            plot(xc1, y1, type = "p", ylim = c(min(y1, y2), max(y1, 
                y2)), xlim = c(min(xc1, xc2), max(xc1, xc2)), 
                log = "y", xlab = "Magnitude", ylab = "Cumulative number of events", 
                pch = 2, cex.lab = 1.3)
            points(xc2, y2, pch = 1)
            ind2 = which(round(xc2, 6) == m02)
            A2 = log10(y2[ind2]) + b2 * m02
        }
        else {
            if (oplt == "p") {
                plot(xc1, y1, type = "p", ylim = c(min(y1), max(y1)), 
                  xlim = c(min(xc1), max(xc1)), log = "y", xlab = "Magnitude", 
                  ylab = "Cumulative number of events", pch = 2, 
                  cex.lab = 1.3)
            }
            else {
                plot(xc1, y1, type = "n", ylim = c(min(y1), max(y1)), 
                  xlim = c(min(xc1), max(xc1)), log = "y", xlab = "Magnitude", 
                  ylab = "Cumulative number of events", pch = 2, 
                  cex.lab = 1.3)
                bwfanc = 0.05
                xcmoins = xc1 - bwfanc
                xcplus = xc1 + bwfanc
                xfanc = sort(c(xcmoins, xcplus))
                yfanc = rep(y1, each = 2)
                lines(xfanc, yfanc, lwd = 2)
            }
        }
        segments(m01, 10^(A1 - b1 * m01), max(xc1), 10^(A1 - 
            b1 * max(xc1)), lwd = 3)
        if (!is.null(finame2)) {
            segments(m02, 10^(A2 - b2 * m02), max(xc2), 10^(A2 - 
                b2 * max(xc2)), lwd = 2)
            dm = max(xc1, xc2) - min(xc1, xc2)
            dn = max(y1, y2) - min(y1, y2)
            posx = min(xc1, xc2) + 0.44 * dm
            posy = min(y1, y2) + 0.75 * dn
            text(posx, posy - 0.25 * dn, paste("bootstrap p-value = ", 
                formatC(randomisation$val$p.value, digits = 2, 
                  format = "g")), font = 1, cex = 1.3, pos = 4)
            text(posx, posy - 0.51 * dn, paste("b (circle) = ", 
                formatC(round(b2, 2), digits = 2, width = 4, 
                  format = "f"), "   ", formatC(round(sd2, 2), 
                  digits = 2, width = 4, format = "f")), font = 1, 
                cex = 1.3, pos = 4)
            text(posx + 0.31 * dm, posy - 0.51 * dn, expression(" " %+-% 
                " "), font = 1, cex = 1.3, pos = 4)
        }
        else {
            dm = max(xc1) - min(xc1)
            dn = max(y1) - min(y1)
            posx = min(xc1) + 0.44 * dm
            posy = min(y1) + 0.75 * dn
        }
        text(posx, posy, title, font = 2, cex = 2, pos = 4)
        if (!is.null(finame2)) 
            text(posx, posy - 0.4 * dn, paste("b (triangle) = ", 
                formatC(round(b1, 2), digits = 2, width = 4, 
                  format = "f"), "   ", formatC(round(sd1, 2), 
                  digits = 2, width = 4, format = "f")), font = 1, 
                cex = 1.3, pos = 4)
        if (is.null(finame2)) 
            text(posx, posy - 0.4 * dn, paste("  b - value  = ", 
                formatC(round(b1, 2), digits = 2, width = 4, 
                  format = "f"), "   ", formatC(round(sd1, 2), 
                  digits = 2, width = 4, format = "f")), font = 1, 
                cex = 1.3, pos = 4)
        text(posx + 0.35 * dm, posy - 0.4 * dn, expression(" " %+-% 
            " "), font = 1, cex = 1.3, pos = 4)
        dev.off()
        cat("Plot in", figname, "\n\n")
    }
}
