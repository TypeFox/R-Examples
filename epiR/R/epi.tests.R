"epi.tests" <- function(dat, conf.level = 0.95) {

    ## Create a list to hold all variables:
    elements <- list()

    ## Do all data manipulation within the list:
    elements <- within(elements, {
        N. <- 1 - ((1 - conf.level) / 2)
        z <- qnorm(N., mean = 0, sd = 1)

        ## Exact binomial confidence limits from D. Collett (1999) Modelling binary data. Chapman & Hall/CRC, Boca Raton Florida, p. 24.
        .funincrisk <- function(cdat, conf.level){
            N. <- 1 - ((1 - conf.level) / 2)
            a <- cdat[,1]
            n <- cdat[,2]
            b <- n - a
            p <- a / n

            a. <- ifelse(a == 0, a + 1, a); b. <- ifelse(b == 0, b + 1, b)
            low <- a. /(a. + (b. + 1) * (1 / qf(1 - N., 2 * a., 2 * b. + 2)))
            up <- (a. + 1) / (a. + 1 + b. / (1 / qf(1 - N., 2 * b., 2 * a. + 2)))
            low <- ifelse(a == 0, 0, low)
            up <- ifelse(a == n, 1, up)
            rval <- data.frame(est = p, lower = low, upper = up)
            rval
        }

        ## From Greg Snow, R-sig-Epi, 3 Mar 2008:
        ## My prefered approach (not the only one), is to use the Bayesian interval using a uniform prior (beta(1,1) distribution)
        ## with the binomial (it is easier to do than it looks). Basically find the HPD interval from a beta distribution with parameters s+1 and f+1,
        ## where s and f are successes (correct test results) and failures (incorrect test results).

        ## I use the hpd function from the TeachingDemos package, but there are others as well (I'm a bit biased towards that package).

        ## For example, to calculate the 95% confidence interval for sensitivity when you have 95 true positives and 5 false negatives you would just
        ## type (after installing and loading the package):
        ## hpd(qbeta, shape1 = 96, shape2 = 6)

        ## And the 2 numbers are limits of a 95% confidence interval. I like this approach because it still gives sensible results when you
        ## have no false negatives (or false positives for specificity).

        ##  hpd. <- function(posterior.icdf, conf = conf.level, tol = 1e-08, ...){
        ##  conf <- min(conf, 1 - conf)
        ##  f <- function(x, posterior.icdf, conf, ...) {
        ##  posterior.icdf(1 - conf + x, ...) - posterior.icdf(x, ...)
        ##  }
        ##  out <- optimize(f, c(0, conf), posterior.icdf = posterior.icdf, conf = conf, tol = tol, ...)
        ##  return(c(posterior.icdf(out$minimum, ...), posterior.icdf(1 - conf + out$minimum, ...)))
        ## }

        ## =================
        ## DECLARE VARIABLES
        ## =================

        ## --------| D+ --| D- --| Total
        ## Test +  | a    | b    | N1
        ## Test -  | c    | d    | N0
        ## --------|------|------|------
        ## Total   | M1   | M0   | total

        a <- dat[1]
        b <- dat[3]
        c <- dat[2]
        d <- dat[4]

        ## Total disease pos:
        M1 <- a + c
        ## Total disease neg:
        M0 <- b + d
        ## Total test pos:
        N1 <- a + b
        ## Total test neg:
        N0 <- c + d
        ## Total subjects:
        total <- a + b + c + d

        ## True prevalence:
        tdat <- as.matrix(cbind(M1, total))
        trval <- .funincrisk(tdat, conf.level)
        tp <- trval$est; tp.low <- trval$lower; tp.up <- trval$upper

        ## Greg Snow:
        ## r <- M1; n <- total
        ## alpha1 <- r + 1
        ## alpha2 <- n - r + 1
        ## tp <- r/n
        ## tp.low <- hpd.(qbeta, shape1 = alpha1, shape2 = alpha2)[1]
        ## tp.up <- hpd.(qbeta, shape1 = alpha1, shape2 = alpha2)[2]

        ## Altman:
        ## q <- 1 - p
        ## A <- (2 * r) + (z * z)
        ## B <- z * sqrt((z * z) + (4 * r * q))
        ## C <- 2 * (n + (z * z))
        ## tp <- p
        ## tp.low <- (A - B) / C
        ## tp.up <- (A + B) / C

        tprev <- data.frame(est = tp, lower = tp.low, upper = tp.up)

        ## Apparent prevalence:
        tdat <- as.matrix(cbind(N1, total))
        trval <- .funincrisk(tdat, conf.level)
        ap <- trval$est; ap.low <- trval$lower; ap.up <- trval$upper

        ## Greg Snow:
        ## r <- N1; n <- total
        ## alpha1 <- r + 1
        ## alpha2 <- n - r + 1
        ## ap <- r/n
        ## ap.low <- hpd.(qbeta, shape1 = alpha1, shape2 = alpha2)[1]
        ## ap.up <- hpd.(qbeta, shape1 = alpha1, shape2 = alpha2)[2]

        ## Altman:
        ## q <- 1 - p
        ## A <- (2 * r) + (z * z)
        ## B <- z * sqrt((z * z) + (4 * r * q))
        ## C <- 2 * (n + (z * z))
        ## ap <- p
        ## ap.low <- (A - B) / C
        ## ap.up <- (A + B) / C

        aprev <- data.frame(est = ap, lower = ap.low, upper = ap.up)

        ## Sensitivity:
        tdat <- as.matrix(cbind(a, M1))
        trval <- .funincrisk(tdat, conf.level)
        se <- trval$est; se.low <- trval$lower; se.up <- trval$upper

        ## Greg Snow:
        ## r <- a; n <- M1
        ## alpha1 <- r + 1
        ## alpha2 <- n - r + 1
        ## se <- r/n
        ## se.low <- hpd.(qbeta, shape1 = alpha1, shape2 = alpha2)[1]
        ## se.up <- hpd.(qbeta, shape1 = alpha1, shape2 = alpha2)[2]

        ## Altman:
        ## q <- 1 - p
        ## A <- (2 * r) + (z * z)
        ## B <- z * sqrt((z * z) + (4 * r * q))
        ## C <- 2 * (n + (z * z))
        ## se <- p
        ## se.low <- (A - B) / C
        ## se.up <- (A + B) / C

        sensitivity <- data.frame(est = se, lower = se.low, upper = se.up)

        ## Specificity:
        tdat <- as.matrix(cbind(d, M0))
        trval <- .funincrisk(tdat, conf.level)
        sp <- trval$est; sp.low <- trval$lower; sp.up <- trval$upper

        ## Greg Snow:
        ## r <- d; n <- M0
        ## alpha1 <- r + 1
        ## alpha2 <- n - r + 1
        ## sp <- r/n
        ## sp.low <- hpd.(qbeta, shape1 = alpha1, shape2 = alpha2)[1]
        ## sp.up <- hpd.(qbeta, shape1 = alpha1, shape2 = alpha2)[2]

        ## Altman:
        ## q <- 1 - p
        ## A <- (2 * r) + (z * z)
        ## B <- z * sqrt((z * z) + (4 * r * q))
        ## C <- 2 * (n + (z * z))
        ## sp <- p
        ## sp.low <- (A - B) / C
        ## sp.up <- (A + B) / C

        specificity <- data.frame(est = sp, lower = sp.low, upper = sp.up)

        ## Positive predictive value:
        tdat <- as.matrix(cbind(a, N1))
        trval <- .funincrisk(tdat, conf.level)
        ppv <- trval$est; ppv.low <- trval$lower; ppv.up <- trval$upper

        ## Greg Snow:
        ## r <- a; n <- N1
        ## alpha1 <- r + 1
        ## alpha2 <- n - r + 1
        ## ppv <- r/n
        ## ppv.low <- hpd.(qbeta, shape1 = alpha1, shape2 = alpha2)[1]
        ## ppv.up <- hpd.(qbeta, shape1 = alpha1, shape2 = alpha2)[2]

        ## Altman:
        ## q <- 1 - p
        ## A <- (2 * r) + (z * z)
        ## B <- z * sqrt((z * z) + (4 * r * q))
        ## C <- 2 * (n + (z * z))
        ## ppv <- p
        ## ppv.low <- (A - B) / C
        ## ppv.up <- (A + B) / C

        pv.positive <- data.frame(est = ppv, lower = ppv.low, upper = ppv.up)

        ## Negative predictive value:
        tdat <- as.matrix(cbind(d, N0))
        trval <- .funincrisk(tdat, conf.level)
        npv <- trval$est; npv.low <- trval$lower; npv.up <- trval$upper

        ## Greg Snow:
        ## r <- d; n <- N0
        ## alpha1 <- r + 1
        ## alpha2 <- n - r + 1
        ## npv <- r/n
        ## npv.low <- hpd.(qbeta, shape1 = alpha1, shape2 = alpha2)[1]
        ## npv.up <- hpd.(qbeta, shape1 = alpha1, shape2 = alpha2)[2]

        ## Altman:
        ## q <- 1 - p
        ## A <- (2 * r) + (z * z)
        ## B <- z * sqrt((z * z) + (4 * r * q))
        ## C <- 2 * (n + (z * z))
        ## npv <- p
        ## npv.low <- (A - B) / C
        ## npv.up <- (A + B) / C

        pv.negative <- data.frame(est = npv, lower = npv.low, upper = npv.up)

        ## Likelihood ratio of a positive test. Confidence intervals from Simel et al. (1991)
        ## lrpos <- se / (1 - sp)
        lrpos <- (a/M1) / (1 - (d/M0))
        lrpos.low <- exp(log(lrpos) - z * sqrt((1 - se) / (M1 * se) + (sp)/(M0 * (1 - sp))))
        lrpos.up <-  exp(log(lrpos) + z * sqrt((1 - se) / (M1 * se) + (sp)/(M0 * (1 - sp))))

        lr.positive <- data.frame(est = lrpos, lower = lrpos.low, upper = lrpos.up)


        ## Likelihood ratio of a negative test. Confidence intervals from Simel et al. (1991)
        ## lrpos <- se / (1 - sp)
        lrneg <- (1 - (a/M1)) / (d/M0)
        lrneg.low <- exp(log(lrneg) - z * sqrt((se)/(M1 * (1 - se)) + (1 - sp)/(M0 * (sp))))
        lrneg.up <-  exp(log(lrneg) + z * sqrt((se)/(M1 * (1 - se)) + (1 - sp)/(M0 * (sp))))

        lr.negative <- data.frame(est = lrneg, lower = lrneg.low, upper = lrneg.up)


        ## Diagnostic accuracy (from Scott et al. (2008)):
        tdat <- as.matrix(cbind((a + d), total))
        trval <- .funincrisk(tdat, conf.level)
        da <- trval$est; da.low <- trval$lower; da.up <- trval$upper

        ## Greg Snow:
        ## r <- (a + d); n <- total
        ## p <- r/n
        ## alpha1 <- r + 1
        ## alpha2 <- n - r + 1
        ## da <- r/n
        ## da.low <- hpd.(qbeta, shape1 = alpha1, shape2 = alpha2)[1]
        ## da.up <- hpd.(qbeta, shape1 = alpha1, shape2 = alpha2)[2]

        ## Altman:
        ## q <- 1 - p
        ## A <- (2 * r) + (z * z)
        ## B <- z * sqrt((z * z) + (4 * r * q))
        ## C <- 2 * (n + (z * z))
        ## da <- p
        ## da.low <- (A - B) / C
        ## da.up <- (A + B) / C

        diag.acc <- data.frame(est = da, lower = da.low, upper = da.up)


        ## Diagnostic odds ratio (from Scott et al. (2008)):
        dOR.p <- (a * d) / (b * c)
        lndOR <- log(dOR.p)
        lndOR.var <- 1/a + 1/b + 1/c + 1/d
        lndOR.se <- sqrt(1/a + 1/b + 1/c + 1/d)
        lndOR.l <- lndOR - (z * lndOR.se)
        lndOR.u <- lndOR + (z * lndOR.se)
        dOR.se <- exp(lndOR.se)
        dOR.low <- exp(lndOR.l)
        dOR.up <- exp(lndOR.u)

        diag.or <- data.frame(est = dOR.p, lower = dOR.low, upper = dOR.up)

        ## Number needed to diagnose (from Scott et al. (2008)):
        ndx <- 1 / (se - (1 - sp))
        ndx.1 <- 1 / (se.low - (1 - sp.low))
        ndx.2 <- 1 / (se.up - (1 - sp.up))
        ndx.low <- min(ndx.1, ndx.2)
        ndx.up <- max(ndx.1, ndx.2)

        nnd <- data.frame(est = ndx, lower = ndx.low, upper = ndx.up)

        ## Youden's index (from Bangdiwala et al. (2008)):
        c.p <- se - (1 - sp)
        c.1 <- se.low - (1 - sp.low)
        c.2 <- se.up - (1 - sp.up)
        c.low <- min(c.1, c.2)
        c.up <- max(c.1, c.2)

        youden <- data.frame(est = c.p, lower = c.low, upper = c.up)

    })

    rval <- list(
        aprev    = elements$aprev,
        tprev    = elements$tprev,
        se       = elements$sensitivity,
        sp       = elements$specificity,
        diag.acc = elements$diag.acc,
        diag.or  = elements$diag.or,
        nnd      = elements$nnd,
        youden   = elements$youden,
        ppv      = elements$pv.positive,
        npv      = elements$pv.negative,
        plr      = elements$lr.positive,
        nlr      = elements$lr.negative)

    ## Define tab:
    r1 <- with(elements, c(a, b, N1))
    r2 <- with(elements, c(c, d, N0))
    r3 <- with(elements, c(M1, M0, M0 + M1))
    tab <- as.data.frame(rbind(r1, r2, r3))
    colnames(tab) <- c("   Outcome +", "   Outcome -", "     Total")
    rownames(tab) <- c("Test +", "Test -", "Total")
    tab <- format.data.frame(tab, digits = 3, justify = "right")

    out <- list(conf.level = conf.level, elements = elements, rval = rval, tab = tab)

    class(out) <- "epi.tests"
    return(out)
}


## Print method for epi.tests:
print.epi.tests <- function(x, ...) {

    print(x$tab, ...)
    cat("\nPoint estimates and", x$conf.level * 100, "%", "CIs:")
    cat("\n---------------------------------------------------------")

    with(x$rval, {

        cat(sprintf("\nApparent prevalence                    %.2f (%.2f, %.2f)",
                    aprev$est,
                    aprev$lower,
                    aprev$upper
                    ))
        cat(sprintf("\nTrue prevalence                        %.2f (%.2f, %.2f)",
                    tprev$est,
                    tprev$lower,
                    tprev$upper
                    ))

        cat(sprintf("\nSensitivity                            %.2f (%.2f, %.2f)",
                    se$est,
                    se$lower,
                    se$upper
                    ))

        cat(sprintf("\nSpecificity                            %.2f (%.2f, %.2f)",
                    sp$est,
                    sp$lower,
                    sp$upper
                    ))

        cat(sprintf("\nPositive predictive value              %.2f (%.2f, %.2f)",
                    ppv$est,
                    ppv$lower,
                    ppv$upper
                    ))

        cat(sprintf("\nNegative predictive value              %.2f (%.2f, %.2f)",
                    npv$est,
                    npv$lower,
                    npv$upper
                    ))

        cat(sprintf("\nPositive likelihood ratio              %.2f (%.2f, %.2f)",
                    plr$est,
                    plr$lower,
                    plr$upper
        ))
        
        cat(sprintf("\nNegative likelihood ratio              %.2f (%.2f, %.2f)",
                    nlr$est,
                    nlr$lower,
                    nlr$upper
        ))
        
    })
    cat("\n---------------------------------------------------------")
    cat("\n")
}


## Summary method for epi.tests:

summary.epi.tests <- function(object, ...) {

    ## Create a data frame:
    out <- do.call(rbind, object$rval)
    
    ## Return it:
    return(out)
}