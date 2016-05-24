gx.md.gait <-
function (xx, wts = NULL, trim = -1, mvtstart = FALSE, mcdstart = FALSE, 
    main = deparse(substitute(xx)), ifadd = c(0.98, 0.95, 0.9), 
    cexf = 0.6, cex = 0.8, ...) 
{
    if (!is.matrix(xx)) 
        stop("  ", deparse(substitute(xx)), " is not a Matrix")
    temp.x <- remove.na(xx)
    x <- temp.x$x
    n <- temp.x$n
    p <- temp.x$m
    matnames <- dimnames(xx)
    matnames[[1]] <- c(1:n)
    if (mcdstart) {
        save <- cov.mcd(x)
        wts <- rep(0, n)
        wts[save$best] <- 1
        xmean <- save$center
        xcov <- save$cov
        proc <- "mcd"
    }
    else {
        if (is.null(wts)) 
            wts <- rep(1, n)
        save <- cov.wt(x, wt = wts)
        xmean <- save$center
        xcov <- save$cov
        proc <- " "
    }
    md <- mahalanobis(x, xmean, xcov)
    frame()
    if (trim > 0) {
        oldpar <- par()
        on.exit(par(oldpar))
        par(mfrow = c(1, 2))
        if (proc == "mcd") 
            wts <- rep(1, n)
        ncore.old <- sum(wts)
        if (trim < 1) {
            ntrim <- floor(ncore.old * trim)
            ptrim <- trim
            proc <- "mvt"
        }
        else {
            if (mvtstart) {
                ncore.old <- n
                wts <- rep(1, n)
            }
            ntrim <- trim
            ptrim <- -1
            proc <- "gait"
        }
        ncore.new <- ncore.old - ntrim
        if (ncore.new < 5 * p) 
            cat("  *** Proceed with Care, ncore < 5p ***\n")
        if (ncore.new < 3 * p) 
            cat("  *** Proceed with Great Care, ncore < 3p ***\n")
        n1 <- ncore.new + 1
        ordered <- order(md)
        for (i in n1:n) wts[ordered[i]] <- 0
        save <- cov.wt(x, wt = wts)
        xmean <- save$center
        xcov <- save$cov
        xsd <- sqrt(diag(xcov))
        md <- mahalanobis(x, xmean, xcov)
        new.core.md <- md[wts[1:n] == 1]
        gx.md.plt0(new.core.md, ncore.new, p, trim = 0, ptrim = ptrim, 
            proc = "", main = "Core Subset", ifadd = ifadd, cexf = cexf, 
            cex = cex, ...)
        gx.md.plt0(md, n, p, trim = n - ncore.new, ptrim = ptrim, 
            proc = proc, main = main, ifadd = ifadd, cexf = cexf, 
            cex = cex, ...)
        par(mfrow = c(1, 1))
    }
    else {
        if (n <= 5 * p) 
            cat("  *** Proceed with Care, n < 5p ***\n")
        gx.md.plt0(md, n, p, trim = 0, proc = proc, main = main, 
            ifadd = ifadd, cexf = cexf, cex = cex, ...)
        ncore.new <- sum(wts)
        ptrim <- -1
        xsd <- sqrt(diag(xcov))
    }
    temp <- (ncore.new - p)/(p * (ncore.new + 1))
    ppm <- 1 - pf(temp * md, p, ncore.new - p)
    invisible(list(main = main, input = deparse(substitute(xx)), 
        matnames = matnames, proc = proc, wts = wts, n = n, nc = ncore.new, 
        p = p, ifilr = FALSE, ptrim = ptrim, mean = xmean, cov = xcov, 
        sd = xsd, md = md, ppm = ppm))
}
