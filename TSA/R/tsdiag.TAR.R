tsdiag.TAR <-
function (object, gof.lag,  col = "red", 
    xlab = "t", ...) 
{
    box.ljung <- function(x, xy1, xy2, nlag = 1, is.print = TRUE) {
        x = x/var(x)^0.5
        em = NULL
        for (i in 1:nlag) {
            em = cbind(em, zlag(x, i))
        }
        em[is.na(em)] = 0
        n = dim(xy1)[1]
        H = t(em) %*% cbind(xy1, xy2)
        dimH = dim(H)
        denom = matrix(rep(n - seq(nlag), dimH[2]), nrow = dimH[1], 
            ncol = dimH[2])
        H = H/denom
        D1 = solve(t(xy1) %*% xy1/n)
        D2 = solve(t(xy2) %*% xy2/n)
        dim1 = dim(D1)
        dim2 = dim(D2)
        D = rbind(cbind(D1, matrix(0, nrow = dim1[1], ncol = dim2[2])), 
            cbind(matrix(0, nrow = dim2[1], ncol = dim1[2]), 
                D2))
        temp = t(em) - H %*% D %*% t(cbind(xy1, xy2))
        Sigma = temp %*% t(temp)/n
        eig = eigen(Sigma)
        eigv = eig$values
        eigm = eig$vectors
        eigv[eigv < 10^(-6) * max(eigv)] = 0
        eigv[eigv > 0] = sqrt(1/eigv[eigv > 0])
        a <- diag(eigv) %*% t(eigm) %*% acf(x, lag.max = nlag, 
            plot = FALSE)$acf
        Q = n * sum(a^2)
        df <- sum(eigv > 0)
        p.value = signif(1 - pchisq(Q, df = df), 4)
        Q1 <- if (is.print) {
            cat1(is.print = is.print, "\n Box-Ljung statistic =", 
                signif(Q, 4), " on ", round(df), " df with\np-value = ", 
                p.value = signif(1 - pchisq(Q, df = df), 4), 
                "\n")
        }
        invisible(list(Q = Q, p.value = p.value, Q1 = (n * sum(acf(x, 
            lag.max = nlag, plot = FALSE)$acf)^2)))
    }
    cat1 <- function(..., is.print = TRUE, file = "", sep = " ", 
        fill = FALSE, labels = NULL, append = FALSE) {
        if (is.print) {
            if (is.character(file)) 
                if (file == "") 
                  file <- stdout()
                else if (substring(file, 1, 1) == "|") {
                  file <- pipe(substring(file, 2), "w")
                  on.exit(close(file))
                }
                else {
                  file <- file(file, ifelse(append, "a", "w"))
                  on.exit(close(file))
                }
            cat(..., file=file, sep=sep, fill=fill, labels=labels, append=append)
        }
        invisible(NULL)
    }
    opar=par(mfrow = c(3, 1), mar = c(3, 4, 3, 2) + 0.1, oma = c(1, 
        0, 2, 0))
    dxy1 = object$dxy1
    dxy2 = object$dxy2
    sort.l = object$sort.l
    dxy1[sort.l, ] = dxy1
    dxy2[sort.l, ] = dxy2
    std.res = object$std.res
    p1 = object$p1
    p2 = object$p2
    d = object$d
    thd = object$thd
    qr1 = object$qr1
    qr2 = object$qr2
    n1 = object$n1
    n2 = object$n2
    n = n1 + n2
    if (missing(gof.lag)) 
        nlag = min(length(std.res)/4, 24)
    else nlag = gof.lag
    plot(std.res, xlab = xlab, ylab = "Standardized Residuals", 
        ...)
    limit = qnorm(0.25/n)
    abline(h = 0)
    abline(h = limit, col = col, lty = 2)
    abline(h = -limit, col = col, lty = 2)
    acf(std.res, xlab = "Lag", ylab = "ACF of Residuals", main='', ci.col = col, 
        ...)
    BL.v = NULL
    for (i in 1:nlag) {
        BL.v = c(BL.v, box.ljung(std.res, dxy1, dxy2, nlag = i, 
            is.print = FALSE)$p.value)
    }
    plot(y = BL.v, x = 1:nlag, ylab = "P-values", xlab = "Number of Lags", 
        ylim = range(c(BL.v, 0)), ...)
    abline(h = 0.05, lty = 2, col = col)
    par(opar)
    invisible()
}
