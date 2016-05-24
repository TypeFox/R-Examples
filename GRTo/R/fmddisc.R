fmddisc <-
function (file, header = FALSE, colid = 1, nrep = 200, title = "fmddisc") 
{
    delta = 0.1
    "mbass" <- function(a, delta = 0.1, plot = TRUE, alldisc = FALSE, 
        bs = 0) {
        mba <- function(x) {
            fmbass(x, delta, plot, alldisc)
        }
        if (bs == 0) {
            res <- mba(a)
        }
        else {
            res = bootstrap::bootstrap(a, abs(bs), mba)
            if (bs > 0) {
                res = quantile(res$thetastar, probs = c(0.05, 
                  0.5, 0.95), na.rm = TRUE)
            }
        }
        invisible(res)
    }
    "fmbass" <- function(a, delta = 0.1, plot = TRUE, alldisc = FALSE) {
        if (plot) {
            par(mfrow = c(1, 1))
        }
        tau <- numeric()
        pva <- numeric()
        minmag <- min(a, na.rm = TRUE)
        maxmag <- max(a, na.rm = TRUE)
        g_r <- hist(a, plot = FALSE, breaks = round(seq((minmag - 
            delta/2), (maxmag + delta/2), delta), 3))
        n <- length(g_r$density)
        xc <- round(seq(minmag, maxmag, delta)[1:(n - 1)], 3)
        log_nc <- log10((1/delta) * (length(a) - cumsum(g_r$counts)[1:(n - 
            1)]) * delta)
        x <- round(seq(minmag, maxmag, delta), 3)
        log_n <- log10((1/delta) * g_r$counts * delta)
        x <- x[is.finite(log_n)]
        log_n <- log_n[is.finite(log_n)]
        sl <- diff(log_n)/diff(x)
        xsl <- x[2:length(x)]
        if (plot) {
            plot(xc, (10^log_nc), type = "p", ylim = c(1, length(a)), 
                log = "y", xlab = "Magnitude", ylab = "Number of events", 
                pch = 1)
            points(x, (10^log_n), pch = 2)
        }
        niter <- 6
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
            if ((n1[1] > 2) && (n1[1] <= (N - 2)) && (wilcox.test(xn1, 
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
        m0 = c(signif(xsl[tau[ip[1]]]), signif(xsl[tau[ip[2]]]))
        if (alldisc) {
            return(print(list(discmag = xsl[tau], p = v_pva, 
                threshold = m0)))
        }
        invisible(m0)
    }
    b = read.table(file, header = header)
    a <- numeric()
    a <- b[, colid]
    a = round(a, 1)
    if (any(a < -3)) 
        stop("Values out of range. Are you sure that input values are magnitude values ?")
    if (any(a > 9.5)) 
        stop("Values out of range. Are you sure that input values are magnitude values ?")
    minmag <- min(a, na.rm = TRUE)
    maxmag <- max(a, na.rm = TRUE)
    figname = paste(title, "_disc.png", sep = "")
    resname1 = paste(title, "_disc1.res", sep = "")
    resname2 = paste(title, "_disc2.res", sep = "")
    cat("\n...you just have to wait...\n\n...Generating", nrep, 
        "replicates...\n\n")
    png(filename = figname, width = 2400, height = 2715, res = 360)
    g_r <- hist(a, plot = FALSE, breaks = round(seq((minmag - 
        delta/2), (maxmag + delta/2), delta), 3))
    n <- length(g_r$density)
    x <- round(seq(minmag, maxmag, delta), 3)
    log_n1 <- log10((1/delta) * (length(a) - cumsum(g_r$counts)) * 
        delta)
    log_n2 <- log10((1/delta) * g_r$counts * delta)
    y = x[is.finite(log_n2)]
    x = x[is.finite(log_n1)]
    log_n1 = log_n1[is.finite(log_n1)]
    log_n2 = log_n2[is.finite(log_n2)]
    par(mar = c(4.1, 4.1, 2.1, 1.1), font = c(2), font.main = c(2), 
        font.axis = c(2), font.lab = c(2))
    layout(matrix(c(1, 2, 3), 3, 1, byrow = FALSE), widths = c(lcm(9.1)), 
        heights = c(lcm(6.1), lcm(5.8), lcm(5.8)), respect = TRUE)
    par(xaxs = "i")
    plot(x, (10^log_n1), type = "p", ylim = c(1, length(a)), 
        log = "y", xlab = "Magnitude", ylab = "Number of events", 
        cex = 1.6, pch = 1, main = title)
    points(y, (10^log_n2), cex = 2.1, pch = 2)
    text(0.75 * (min(x) + max(x)), 0.5 * 10^max(log_n1), "(A) inc. and cum. FMDs")
    par(xaxs = "r")
    b = read.table(file, header = header)
    a <- numeric()
    a <- b[, colid]
    a = round(a, 1)
    if (any(a < -3)) 
        stop("Values out of range. Are you sure that input values are magnitude values ?")
    if (any(a > 9.5)) 
        stop("Values out of range. Are you sure that input values are magnitude values ?")
    y = mbass(a, delta = delta, plot = FALSE, alldisc = FALSE, 
        bs = -nrep)
    cat("\n***********************\nMAIN DISCONTINUITY\n")
    cat("\nQuantiles of the main discontinuity in the frequency-magnitude distribution:\n")
    print(quantile(y$thetastar[1, ], probs = c(0.05, 0.5, 0.95), 
        na.rm = TRUE))
    cat("Number of valid replicates:\n")
    print(sum(!is.na(y$thetastar[1, ])))
    cat("\nBootstrap mean and bootstrap standard-error (standard deviation of bootstrap replicate estimates):\n")
    print(round(mean(y$thetastar[1, ], na.rm = TRUE), 3))
    print(round(sd(y$thetastar[1, ], na.rm = TRUE), 3))
    cat("Bootstrap margin of errors at the 90% normal confidence level:\n")
    print(round(1.645 * sd(y$thetastar[1, ], na.rm = TRUE), 3))
    cat("\n***********************\nAUXILIARY DISCONTINUITY\n")
    cat("\nQuantiles of the auxiliary discontinuity in the frequency-magnitude distribution:\n")
    print(quantile(y$thetastar[2, ], probs = c(0.05, 0.5, 0.95), 
        na.rm = TRUE))
    cat("Number of valid replicates:\n")
    print(sum(!is.na(y$thetastar[2, ])))
    cat("\nBootstrap mean and bootstrap standard-error (standard deviation of bootstrap replicate estimates):\n")
    print(round(mean(y$thetastar[2, ], na.rm = TRUE), 3))
    print(round(sd(y$thetastar[2, ], na.rm = TRUE), 3))
    cat("Bootstrap margin of errors at the 90% normal confidence level:\n")
    print(round(1.645 * sd(y$thetastar[2, ], na.rm = TRUE), 3))
    hist(y$thetastar[1, ], main = "(B) Main discontinuity", xlab = "Discontinuity magnitude", 
        ylim = c(0, nrep), breaks = round(seq((minmag - delta/2), 
            (maxmag + delta/2), delta), 3))
    hist(y$thetastar[2, ], main = "(C) Auxiliary discontinuity", 
        xlab = "Discontinuity magnitude", ylim = c(0, nrep), 
        breaks = round(seq((minmag - delta/2), (maxmag + delta/2), 
            delta), 3))
    dev.off()
    write.table(file = resname1, y$thetastar[1, ], row.names = FALSE, 
        col.names = FALSE)
    write.table(file = resname2, y$thetastar[2, ], row.names = FALSE, 
        col.names = FALSE)
    cat("\nReplicates for magnitude discontinuities in files", 
        resname1, "(main) and\n", resname2, "(auxiliary)\n")
    cat("\nFmddisc plot in", figname, "\n\n")
    disc = list()
    disc$quant1 = quantile(y$thetastar[1, ], probs = c(0.05, 
        0.5, 0.95), na.rm = TRUE)
    disc$valid[1] = sum(!is.na(y$thetastar[1, ]))
    disc$bmean[1] = round(mean(y$thetastar[1, ], na.rm = TRUE), 
        3)
    disc$bse[1] = round(sd(y$thetastar[1, ], na.rm = TRUE), 3)
    disc$bme[1] = round(1.645 * sd(y$thetastar[1, ], na.rm = TRUE), 
        3)
    disc$quant2 = quantile(y$thetastar[2, ], probs = c(0.05, 
        0.5, 0.95), na.rm = TRUE)
    disc$valid[2] = sum(!is.na(y$thetastar[2, ]))
    disc$bmean[2] = round(mean(y$thetastar[2, ], na.rm = TRUE), 
        3)
    disc$bse[2] = round(sd(y$thetastar[2, ], na.rm = TRUE), 3)
    disc$bme[2] = round(1.645 * sd(y$thetastar[2, ], na.rm = TRUE), 
        3)
    disc
}
