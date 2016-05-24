## last modified May 2008

plot.mix <- function (x, mixpar = NULL, dist = "norm", root = FALSE, 
    ytop = NULL, clwd = 1, main, sub, xlab, ylab, bty, BW = FALSE, 
    ...) 
{
    mixobj<-x
    if (inherits(mixobj, "mix")) {
        mixdat <- mixobj$mixdata
        mixpar <- mixobj$parameters
        dist <- mixobj$distribution
    }
    else mixdat <- mixobj
    black <- "black"
    blue <- "blue"
    red <- "red"
    forestgreen <- "forestgreen"
    gray = "gray"
    darkgray = "black"
    lty1 <- 1
    lty2 <- 2
    if (missing(xlab)) 
        xlab <- dimnames(mixdat)[[2]][1]
    if (missing(ylab)) {
        if (root) 
            ylab <- "Square Root of Probability Density"
        else ylab <- "Probability Density"
    }
    if (missing(bty)) 
        bty <- "o"
    ntot <- sum(mixdat[, 2])
    m <- nrow(mixdat)
    iwid <- mixdat[2:(m - 1), 1] - mixdat[1:(m - 2), 1]
    iwid <- c(2 * iwid[1], iwid, 2 * iwid[m - 2])
    if (mixdat[1, 1] > 0 & !is.na(match(dist, c("lnorm", "gamma", 
        "weibull")))) 
        iwid[1] <- min(iwid[1], mixdat[1, 1])
    if (mixdat[1, 1] > 0 & !is.na(match(dist, c("binom", "nbinom", 
        "pois")))) 
        if (mixdat[1, 1] == min(iwid[1], mixdat[1, 1])) 
            iwid[1] <- mixdat[1, 1] + 0.5
    idens <- (mixdat[, 2]/iwid)/ntot
    hx <- rep(c(mixdat[1, 1] - iwid[1], mixdat[-m, 1], mixdat[m - 
        1, 1] + iwid[m]), rep(2, m + 1))
    hy <- c(0, rep(idens, rep(2, m)), 0)
    if (!is.null(mixpar)) {
        k <- nrow(mixpar)
        if (ncol(mixpar) > 3) 
            par3 <- mixpar[, 4]
        else par3 <- rep(0, k)
        prop <- mixpar[, 1]
        if (dist == "norm") {
            par1 <- mixpar[, 2]
            par2 <- mixpar[, 3]
            if (missing(sub)) 
                sub <- "Normal Mixture"
        }
        else if (dist == "lnorm") {
            par2 <- sqrt(log((mixpar[, 3]/(mixpar[, 2] - par3))^2 + 
                1))
            par1 <- log(mixpar[, 2] - par3) - (par2^2)/2
            if (missing(sub)) 
                sub <- "Lognormal Mixture"
        }
        else if (dist == "gamma") {
            par1 <- ((mixpar[, 2] - par3)/mixpar[, 3])^2
            par2 <- (mixpar[, 2] - par3)/(mixpar[, 3]^2)
            if (missing(sub)) 
                sub <- "Gamma Mixture"
        }
        else if (dist == "weibull") {
            par <- weibullpar(mixpar[, 2], mixpar[, 3], par3)
            par1 <- par$shape
            par2 <- par$scale
            if (missing(sub)) 
                sub <- "Weibull Mixture"
        }
        else if (dist == "binom") {
            par1 <- ceiling(mixpar[, 2]^2/(mixpar[, 2] - mixpar[, 
                3]^2))
            par2 <- 1 - mixpar[, 3]^2/mixpar[, 2]
            if (missing(sub)) 
                sub <- "Binomial Mixture"
        }
        else if (dist == "nbinom") {
            par1 <- mixpar[, 2]^2/(mixpar[, 3]^2 - mixpar[, 2])
            if (missing(sub)) 
                sub <- "Negative Binomial Mixture"
        }
        else if (dist == "pois") {
            par1 <- mixpar[, 2]
            if (missing(sub)) 
                sub <- "Poisson Mixture"
        }
        else {
            par1 <- mixpar[, 2]
            par2 <- mixpar[, 3]
            warning(paste("Unknown distribution ", dist, ", using normal", 
                sep = ""))
            dist <- "norm"
        }
        if (root) {
            mwid <- c(mixdat[1, 1] - iwid[1], mixdat[-m, 1]) + 
                iwid/2
            rf <- 2
        }
        else {
            mwid <- 0
            rf <- 1
        }
        lgr <- seq(hx[1], hx[length(hx)], length = 400)
        if (!is.na(match(dist, c("binom", "nbinom", "pois")))) 
            lgr <- floor(lgr[1]):ceiling(lgr[length(lgr)])
        rc <- TRUE
        for (j in 1:rf) {
            pdf <- matrix(0, nrow = k, ncol = length(lgr))
            for (i in 1:k) {
                if (dist == "norm") 
                  pdfi <- dnorm(lgr, par1[i], par2[i])
                else if (dist == "lnorm") 
                  pdfi <- dlnorm(lgr - par3[i], par1[i], par2[i])
                else if (dist == "gamma") 
                  pdfi <- dgamma(lgr - par3[i], par1[i], par2[i])
                else if (dist == "weibull") 
                  pdfi <- dweibull(lgr - par3[i], par1[i], par2[i])
                else if (dist == "binom") 
                  pdfi <- dbinom(round(lgr, 0), par1[i], par2[i])
                else if (dist == "nbinom") 
                  pdfi <- dnbinom(round(lgr, 0), par1[i], mu = mixpar[i, 
                    2])
                else if (dist == "pois") 
                  pdfi <- dpois(round(lgr, 0), par1[i])
                pdf[i, ] <- prop[i] * pdfi
            }
            if (rc) {
                lgr1 <- lgr
                lgr <- mwid
                pdf1 <- pdf
                pdfm1 <- apply(pdf1, 2, sum)
                rc <- FALSE
            }
        }
        if (root) {
            pdfm <- apply(pdf, 2, sum)
            dif <- sqrt(pdfm) - sqrt(idens)
            hy1 <- c(0, rep(pdfm, rep(2, m)), 0)
        }
    }
    if (root) {
        if (is.null(mixpar)) {
            if (is.null(ytop)) 
                ylim <- NULL
            else ylim <- c(0, ytop)
        }
        else {
            if (is.null(ytop)) 
                ylim <- c(min(dif, 0), max(sqrt(hy)))
            else ylim <- c(min(dif, 0), ytop)
        }
        plot(hx, sqrt(hy), type = "n", ylim = ylim, col = black, 
            xlab = xlab, ylab = ylab, bty = bty, ...)
    }
    else {
        if (is.null(ytop)) 
            ylim <- NULL
        else ylim <- c(0, ytop)
        plot(hx, hy, type = "n", ylim = ylim, col = black, xlab = xlab, 
            ylab = ylab, bty = bty, ...)
    }
    title(main = ifelse(missing(main), "", main), sub = ifelse(missing(sub), 
        "", sub))
    if (root & is.null(mixpar)) 
        lines(hx, sqrt(hy), col = ifelse(BW, black, blue))
    else if (root & !is.null(mixpar)) {
        lines(hx, sqrt(hy1), col = ifelse(BW, black, blue))
        rect(c(mixdat[1, 1] - iwid[1], mixdat[-m, 1]), rep(0, 
            m), c(mixdat[-m, 1], mixdat[m - 1, 1] + iwid[m]), 
            dif, col = ifelse(BW, gray, blue))
    }
    else if (!root) 
        lines(hx, hy, col = ifelse(BW, black, blue))
    lines(hx[c(1, length(hx))], c(0, 0), col = black)
    if (!is.null(mixpar)) {
        for (i in 1:k) {
            if (root) 
                lines(lgr1, sqrt(pdf1[i, ]), col = ifelse(BW, 
                  darkgray, red), lty = ifelse(BW, lty2, 1), 
                  lwd = 1 * clwd)
            else lines(lgr1, pdf1[i, ], col = ifelse(BW, darkgray, 
                red), lty = ifelse(BW, lty2, 1), lwd = 1 * clwd)
            points(mixpar[i, 2], 0, pch = 17, col = ifelse(BW, 
                darkgray, red))
        }
        if (root) 
            lines(lgr1, sqrt(pdfm1), col = ifelse(BW, black, 
                forestgreen), lty = ifelse(BW, lty1, 1), lwd = 2 * 
                clwd)
        else lines(lgr1, pdfm1, col = ifelse(BW, black, forestgreen), 
            lty = ifelse(BW, lty1, 1), lwd = 2 * clwd)
    }
}
