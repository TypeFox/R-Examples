predIntPois <-
function (x, k = 1, n.sum = 1, method = "conditional", pi.type = "two-sided", 
    conf.level = 0.95, round.limits = TRUE) 
{
    if (length(k) != 1 || !is.numeric(k) || k != trunc(k) || 
        k < 1 || length(n.sum) != 1 || !is.numeric(n.sum) || 
        n.sum != trunc(n.sum) || n.sum < 1) 
        stop("'k' and 'n.sum' must be a positive integers")
    method <- match.arg(method, c("conditional", "conditional.approx.normal", 
        "conditional.approx.t", "normal.approx"))
    if (method == "conditional" && k > 1) 
        stop("The 'conditonal' method is only implemented for 'k'=1")
    pi.type <- match.arg(pi.type, c("two-sided", "lower", "upper"))
    if (!is.numeric(conf.level) || length(conf.level) > 1 || 
        conf.level <= 0 || conf.level >= 1) 
        stop("'conf.level' must be a scalar greater than 0 and less than 1.")
    if (x.is.est.obj <- data.class(x) == "estimate" || data.class(x) == 
        "estimateCensored") {
        if (x$distribution != "Poisson") 
            stop(paste("'predIntPois' creates prediction intervals", 
                "for a Poisson distribution.  You have supplied an object", 
                "that assumes a different distribution."))
        class.x <- oldClass(x)
        if (!is.null(x$interval)) {
            x <- x[-match("interval", names(x))]
            oldClass(x) <- class.x
        }
        lambda.hat <- x$parameters
        n <- x$sample.size
        ret.list <- x
    }
    else {
        if (!is.vector(x, mode = "numeric")) 
            stop(paste("'x' must be either a list that inherits from", 
                "the class 'estimate', or else a numeric vector"))
        data.name <- deparse(substitute(x))
        if ((bad.obs <- sum(!(x.ok <- is.finite(x)))) > 0) {
            is.not.finite.warning(x)
            x <- x[x.ok]
            warning(paste(bad.obs, "observations with NA/NaN/Inf in 'x' removed."))
        }
        n <- length(x)
        if (n == 0) 
            stop("'x' does not contain any finite, non-missing values")
        if (any(x < 0) || any(x != trunc(x))) 
            stop("All non-missing values of 'x' must be non-negative integers")
        if (all(x == 0)) 
            stop("All finite, non-missing values of 'x' are 0")
        ret.list <- epois(x)
        ret.list$data.name <- data.name
        ret.list$bad.obs <- bad.obs
        lambda.hat <- ret.list$parameters
    }
    sum.x <- n * lambda.hat
    alpha <- 1 - conf.level
    x.hat <- n.sum * lambda.hat
    switch(method, conditional = {
        limits <- ci.normal.approx(theta.hat = x.hat, sd.theta.hat = sqrt((x.hat * 
            (n + n.sum))/n), n = n, df = n - 1, ci.type = pi.type, 
            alpha = alpha, lb = 0)$limits
        lpl <- limits[1]
        upl <- limits[2]
        switch(pi.type, `two-sided` = {
            if (ebinom(0, sum.x, ci = TRUE, ci.type = "two-sided", 
                conf.level = 1 - alpha/2)$interval$limits["UCL"] >= 
                (n.sum/(n.sum + n))) {
                lpl <- 0
                warning(paste("Lower prediction limit not accurate", 
                  "due to discrete nature of Poisson distribution"))
            } else {
                fcn.to.min.two.lpl <- function(lpl, n.weird, 
                  m.weird, y.weird, alpha) (m.weird/(lpl + 1) - 
                  (n.weird/y.weird) * qf(1 - alpha/2, 2 * (lpl + 
                    1), 2 * y.weird))^2
                lpl <- nlminb(start = lpl, objective = fcn.to.min.two.lpl, 
                  lower = 0, n.weird = n, m.weird = n.sum, y.weird = sum.x, 
                  alpha = alpha)$par
                fcn.to.min.two.upl <- function(upl, n.weird, 
                  m.weird, y.weird, alpha) (upl/m.weird - ((y.weird + 
                  1)/n.weird) * qf(1 - alpha/2, 2 * (y.weird + 
                  1), 2 * upl))^2
                upl <- nlminb(start = max(upl, 1), objective = fcn.to.min.two.upl, 
                  lower = 1, n.weird = n, m.weird = n.sum, y.weird = sum.x, 
                  alpha = alpha)$par
            }
        }, lower = {
            if (ebinom(0, sum.x, ci = TRUE, ci.type = "upper", 
                conf.level = conf.level)$interval$limits["UCL"] >= 
                (n.sum/(n.sum + n))) {
                lpl <- 0
                warning(paste("Lower prediction limit not accurate", 
                  "due to discrete nature of Poisson distribution"))
            } else {
                fcn.to.min.one.lpl <- function(lpl, n.weird, 
                  m.weird, y.weird, conf.level) (m.weird/(lpl + 
                  1) - n.weird/y.weird * qf(conf.level, 2 * (lpl + 
                  1), 2 * y.weird))^2
                lpl <- nlminb(start = lpl, objective = fcn.to.min.one.lpl, 
                  lower = 0, n.weird = n, m.weird = n.sum, y.weird = sum.x, 
                  conf.level = conf.level)$par
            }
        }, upper = {
            fcn.to.min.one.upl <- function(upl, n.weird, m.weird, 
                y.weird, conf.level) (upl/m.weird - (y.weird + 
                1)/n.weird * qf(conf.level, 2 * (y.weird + 1), 
                2 * upl))^2
            upl <- nlminb(start = upl, objective = fcn.to.min.one.upl, 
                lower = 0, n.weird = n, m.weird = n.sum, y.weird = sum.x, 
                conf.level = conf.level)$par
        })
        limits <- c(lpl, upl)
    }, conditional.approx.normal = {
        cf <- n.sum/n
        z.crit <- ifelse(pi.type == "two-sided", qnorm(1 - (alpha/k)/2), 
            qnorm(1 - alpha/k))
        hw <- (cf * z.crit^2)/2 + cf * z.crit * sqrt(sum.x * 
            (1 + 1/cf) + z.crit^2/4)
        switch(pi.type, `two-sided` = {
            lpl <- x.hat - hw
            upl <- x.hat + hw
            lpl <- max(0, lpl)
        }, lower = {
            lpl <- x.hat - hw
            lpl <- max(0, lpl)
            upl <- Inf
        }, upper = {
            lpl <- 0
            upl <- x.hat + hw
        })
        limits <- c(lpl, upl)
    }, conditional.approx.t = {
        df <- n - 1
        cf <- n.sum/n
        t.crit <- ifelse(pi.type == "two-sided", qt(1 - (alpha/k)/2, 
            df), qt(1 - alpha/k, df))
        hw <- (cf * t.crit^2)/2 + cf * t.crit * sqrt(sum.x * 
            (1 + 1/cf) + t.crit^2/4)
        switch(pi.type, `two-sided` = {
            lpl <- x.hat - hw
            upl <- x.hat + hw
            lpl <- max(0, lpl)
        }, lower = {
            lpl <- x.hat - hw
            lpl <- max(0, lpl)
            upl <- Inf
        }, upper = {
            lpl <- 0
            upl <- x.hat + hw
        })
        limits <- c(lpl, upl)
    }, normal.approx = {
        if (sum.x <= 10 || x.hat <= 10) warning(paste("Estimated value of 'lambda' and/or", 
            "number of future observations is/are probably too small", 
            "for the normal approximation to work well.\n"))
        limits <- ci.normal.approx(theta.hat = x.hat, sd.theta.hat = sqrt((x.hat * 
            (n + n.sum))/n), n = n, df = n - 1, ci.type = pi.type, 
            alpha = alpha/k, lb = 0)$limits
    })
    if (round.limits) 
        limits <- round(limits, 0)
    names(limits) <- c("LPL", "UPL")
    pi.obj <- list(name = "Prediction", limits = limits, type = pi.type, 
        method = method, conf.level = conf.level, sample.size = n, 
        k = k, n.sum = n.sum)
    oldClass(pi.obj) <- "intervalEstimate"
    ret.list <- c(ret.list, list(interval = pi.obj))
    if (x.is.est.obj) 
        oldClass(ret.list) <- class.x
    else oldClass(ret.list) <- "estimate"
    ret.list
}
