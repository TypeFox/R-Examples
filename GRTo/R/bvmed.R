bvmed <-
function (file, lis, hd = FALSE, colid = 1, nrep = 200, tm = NULL, 
    findtm = TRUE, title = "bvmed") 
{
    delta = 0.1
    "mblm" <- function(formula, dataframe, repeated = TRUE) {
        if (missing(dataframe)) 
            dataframe <- environment(formula)
        term <- as.character(attr(terms(formula), "variables")[-1])
        x = dataframe[[term[2]]]
        y = dataframe[[term[1]]]
        if (length(term) > 2) {
            stop("Only linear models are accepted")
        }
        xx = sort(x)
        yy = y[order(x)]
        n = length(xx)
        slopes = c()
        intercepts = c()
        smedians = c()
        imedians = c()
        if (repeated) {
            for (i in 1:n) {
                slopes = c()
                intercepts = c()
                for (j in 1:n) {
                  if (xx[j] != xx[i]) {
                    slopes = c(slopes, (yy[j] - yy[i])/(xx[j] - 
                      xx[i]))
                    intercepts = c(intercepts, (xx[j] * yy[i] - 
                      xx[i] * yy[j])/(xx[j] - xx[i]))
                  }
                }
                smedians = c(smedians, median(slopes))
                imedians = c(imedians, median(intercepts))
            }
            slope = median(smedians)
            intercept = median(imedians)
        }
        else {
            for (i in 1:(n - 1)) {
                for (j in i:n) {
                  if (xx[j] != xx[i]) {
                    slopes = c(slopes, (yy[j] - yy[i])/(xx[j] - 
                      xx[i]))
                  }
                }
            }
            slope = median(slopes)
            intercepts = yy - slope * xx
            intercept = median(intercepts)
        }
        res = list()
        res$coefficients = c(intercept, slope)
        names(res$coefficients) = c("(Intercept)", term[2])
        res$residuals = y - slope * x - intercept
        names(res$residuals) = as.character(1:length(res$residuals))
        res$fitted.values = x * slope + intercept
        names(res$fitted.values) = as.character(1:length(res$fitted.values))
        if (repeated) {
            res$slopes = smedians
            res$intercepts = imedians
        }
        else {
            res$slopes = slopes
            res$intercepts = intercepts
        }
        res$df.residual = n - 2
        res$rank = 2
        res$terms = terms(formula)
        res$call = match.call()
        res$model = data.frame(y, x)
        res$assign = c(0, 1)
        if (missing(dataframe)) {
            res$effects = lm(formula)$effects
            res$qr = lm(formula)$qr
        }
        else {
            res$effects = lm(formula, dataframe)$effects
            res$qr = lm(formula, dataframe)$qr
        }
        res$effects[2] = sqrt(sum((res$fitted - mean(res$fitted))^2))
        res$xlevels = list()
        names(res$model) = term
        attr(res$model, "terms") = terms(formula)
        class(res) = c("mblm", "lm")
        res
    }
    "fmbass2" <- function(a, delta = 0.1, alldisc = FALSE, tm = NULL, 
        findtm = TRUE) {
        tau <- numeric()
        pva <- numeric()
        minmag <- min(a, na.rm = TRUE)
        g_r <- hist(a, plot = FALSE, breaks = seq((minmag - delta/2), 
            (max(a, na.rm = TRUE) + delta/2), delta))
        n <- length(g_r$density)
        xc <- seq(minmag, max(a, na.rm = TRUE), delta)[1:(n - 
            1)]
        y = numeric()
        for (i in 1:length(xc)) y[i] = sum(a >= xc[i])
        log_nc = log10(y)
        x <- seq(minmag, max(a, na.rm = TRUE), delta)
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
        if (alldisc) 
            return(print(list(discmag = xsl[tau], p = v_pva, 
                threshold = m0)))
        invisible(m0)
        y = log_n
        yc = log_nc
        if (!is.null(tm)) 
            m0[1] = tm
        if (!is.na(m0[1])) {
            x2 = x[which(x >= m0[1])]
            y2 = y[which(x >= m0[1])]
            xc2 = xc[which(xc >= m0[1])]
            yc2 = yc[which(xc >= m0[1])]
            fit = mblm(y2 ~ x2, repeated = TRUE)
            fitcoef = -fit$coefficients[2]
            xval = xc2
            yval = yc2
        }
        if (is.null(tm) & !findtm) {
            fit = mblm(y ~ x, repeated = TRUE)
            fitcoef = -fit$coefficients[2]
            xval = xc
            yval = yc
        }
        if (is.null(tm) & findtm & is.na(m0[1])) {
            fitcoef = NA
            xval = NA
            yval = NA
        }
        return(list(thr = m0, coef = fitcoef, X = xval, Y = yval))
    }
    "mbass2" <- function(a, delta = 0.1, alldisc = FALSE, bs = 0, 
        tm = NULL, findtm = TRUE) {
        mba <- function(x) {
            fmbass2(x, delta, alldisc, tm, findtm)
        }
        res <- if (bs == 0) {
            mba(a)
        }
        else {
            bootstrap::bootstrap(a, bs, mba)
        }
        return(res)
    }
    cat("\nThe R bootstrap package is necessary\n")
    no.file <- missing(file)
    no.lis <- missing(lis)
    a <- if (no.file) {
        as.vector(lis)
    }
    else if (no.lis) {
        x = read.table(file, header = hd)
        x[, colid]
    }
    else {
        stop("must specify either 'file' or 'lis'")
    }
    a = round(a, 1)
    cat("\n...You just have to wait...\n")
    z = mbass2(a, delta = delta, alldisc = FALSE, bs = nrep, 
        tm = tm, findtm = findtm)
    cat("\nSize of data set : ", length(a), "events\n")
    cat("\nM0 ('completeness') value or quantiles :\n")
    mthresh = numeric()
    bvalue = numeric()
    for (i in 1:nrep) {
        mthresh[i] = z$thetastar[[i]]$thr[[1]]
        bvalue[i] = z$thetastar[[i]]$coe
    }
    if (nrep != 0) {
        print(quantile(mthresh, probs = c(0.05, 0.5, 0.95), na.rm = TRUE))
        cat("\nREPLICATES SHOWING UNDETERMINED M0:", sum(is.na(mthresh)), 
            "\n")
        cat("\nb-value or b-value quantiles :\n")
        print(signif(quantile(bvalue, probs = c(0.05, 0.5, 0.95), 
            na.rm = TRUE), digits = 3))
        cat("\n")
    }
    bvm = list()
    bvm$quantm = quantile(mthresh, probs = c(0.05, 0.5, 0.95), 
        na.rm = TRUE)
    bvm$mmed = round(median(mthresh, na.rm = TRUE), 3)
    bvm$quantb = quantile(bvalue, probs = c(0.05, 0.5, 0.95), 
        na.rm = TRUE)
    bvm$valid = sum(!is.na(mthresh))
    z0 = mbass2(a, delta = delta, alldisc = FALSE, bs = 0, tm = tm, 
        findtm = findtm)
    bvm$brm = round(z0$coef, 3)
    bvm$bse = round(sd(bvalue, na.rm = TRUE), 3)
    bvm$bme = round(1.645 * sd(bvalue, na.rm = TRUE), 3)
    z1 = mbass2(a, delta = delta, alldisc = FALSE, bs = 0, tm = NULL, 
        findtm = FALSE)
    figname = paste(title, "_bvmed.png", sep = "")
    png(filename = figname, width = 2400, height = 2715, res = 360)
    mval = z1$X
    nval = z1$Y
    plot(mval, (10^nval), type = "p", ylim = c(1, length(a)), 
        log = "y", xlab = "Magnitude", ylab = "Number of events", 
        pch = 1, cex.lab = 1.5, cex.axis = 1.2)
    b = bvm$brm
    sd = bvm$bse
    print(mval)
    print(10^nval)
    epsi = 1e-05
    nthr = z1$Y[which(abs(z1$X - bvm$mmed) < epsi)]
    A = nthr + b * bvm$mmed
    nmag = sum(a >= bvm$mmed)
    w = bvm$mmed
    segments(bvm$mmed, 10^(A - b * bvm$mmed), max(z1$X), 10^(A - 
        b * max(z1$X)), col = "red", lwd = 2)
    dm = max(mval) - min(mval)
    dn = max(10^nval) - min(10^nval)
    posx = min(mval) + 0.44 * dm
    posy = min(nval) + 0.75 * dn
    text(posx, posy, title, font = 2, cex = 2, pos = 4)
    text(posx, posy - 0.35 * dn, paste("b =", formatC(round(b, 
        2), digits = 2, width = 4, format = "f"), "   ", formatC(round(sd, 
        2), digits = 2, width = 4, format = "f")), font = 1, 
        cex = 1.3, pos = 4)
    text(posx + 0.155 * dm, posy - 0.35 * dn, expression(" " %+-% 
        " "), font = 1, cex = 1.3, pos = 4)
    dev.off()
    cat("\nb-value=", b, "with standard-error", sd, "\nCalculation over", 
        nmag, "magnitude values\n\n")
    cat("\nPlot in", figname, "\n\n")
    return(invisible(bvm))
}
