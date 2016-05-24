ci.qnorm <-
function (p, muhat, sdhat, n, method = c("exact", "normal.approx"), 
    ci.type, alpha, digits = 0) 
{
    method <- match.arg(method)
    if (p == 0.5) {
        method <- "Exact"
        ret.obj <- ci.norm(muhat, sdhat, n, ci.type, alpha)
        ret.obj$parameter <- "Median"
    }
    else {
        df <- n - 1
        conf.level <- 1 - alpha
        if (method == "exact") {
            method <- "Exact"
            ncp <- qnorm(p) * sqrt(n)
            sn <- sqrt(n)
            switch(ci.type, `two-sided` = {
                lcl <- muhat + (qT(p = alpha/2, df = df, ncp = ncp)/sn) * 
                  sdhat
                ucl <- muhat + (qT(p = 1 - alpha/2, df = df, 
                  ncp = ncp)/sn) * sdhat
            }, lower = {
                lcl <- muhat + (qT(p = alpha, df = df, ncp = ncp)/sn) * 
                  sdhat
                ucl <- Inf
            }, upper = {
                lcl <- -Inf
                ucl <- muhat + (qT(p = conf.level, df = df, ncp = ncp)/sn) * 
                  sdhat
            })
            ci.limits <- c(lcl, ucl)
            names(ci.limits) <- c("LCL", "UCL")
        }
        else {
            method <- "Normal Approx"
            zp <- qnorm(p)
            xp <- muhat + zp * sdhat
            sd.xp <- (sdhat/sqrt(n)) * sqrt(1 + zp^2/2)
            ci.limits <- ci.normal.approx(xp, sd.xp, n = n, df = df, 
                ci.type = ci.type, alpha = alpha)$limits
        }
        pct <- round(100 * p, digits)
        ret.obj <- list(name = "Confidence", parameter = paste(pct, 
            number.suffix(pct), " %ile", sep = ""), limits = ci.limits, 
            type = ci.type, method = method, conf.level = conf.level, 
            sample.size = n, dof = df)
    }
    oldClass(ret.obj) <- "intervalEstimate"
    ret.obj
}
