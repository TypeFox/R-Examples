twoSampleLinearRankTestCensored <-
function (x, x.censored, y, y.censored, censoring.side = "left", 
    location.shift.null = 0, scale.shift.null = 1, alternative = "two.sided", 
    test = "logrank", variance = "hypergeometric", surv.est = "prentice", 
    shift.type = "location") 
{
    censoring.side <- match.arg(censoring.side, c("left", "right"))
    alternative <- match.arg(alternative, c("two.sided", "less", 
        "greater"))
    test <- match.arg(test, c("logrank", "tarone-ware", "gehan", 
        "peto-peto", "normal.scores.2", "normal.scores.1", "generalized.sign"))
    variance <- match.arg(variance, c("hypergeometric", "permutation", 
        "asymptotic"))
    surv.est <- match.arg(surv.est, c("prentice", "kaplan-meier", 
        "peto-peto", "altshuler"))
    if (test == "logrank") 
        surv.est <- "altshuler"
    shift.type <- match.arg(shift.type, c("location", "scale"))
    if (variance == "asymptotic" && (test != "peto-peto" || surv.est != 
        "prentice")) 
        stop(paste("The test based on the asymptotic variance is only available for", 
            "the 'peto-peto' test with the 'prentice' estimator of survival,", 
            "i.e., test='peto-peto' and surv.est='prentice'"))
    data.name <- c(deparse(substitute(x)), deparse(substitute(y)))
    names(data.name) <- c("x", "y")
    censoring.name <- c(deparse(substitute(x.censored)), deparse(substitute(y.censored)))
    names(censoring.name) <- c("x", "y")
    if (!is.vector(x, mode = "numeric") || is.factor(x)) 
        stop("'x' must be a numeric vector")
    if (!((is.vector(x.censored, mode = "numeric") && !is.factor(x.censored)) || 
        is.vector(x.censored, mode = "logical"))) 
        stop("'x.censored' must be a logical or numeric vector")
    if (length(x.censored) != length(x)) 
        stop("'x.censored' must be the same length as 'x'")
    if ((bad.obs.x <- sum(!(ok <- is.finite(x) & is.finite(as.numeric(x.censored))))) > 
        0) {
        x <- x[ok]
        x.censored <- x.censored[ok]
        warning(paste(bad.obs.x, "observations with NA/NaN/Inf in 'x' and/or 'x.censored' removed."))
    }
    if (is.numeric(x.censored)) {
        if (!all(x.censored == 0 | x.censored == 1)) 
            stop(paste("When 'x.censored' is a numeric vector, all values of", 
                "'x.censored' must be 0 (not censored) or 1 (censored)."))
        x.censored <- as.logical(x.censored)
    }
    n.x.cen <- sum(x.censored)
    if (!is.vector(y, mode = "numeric") || is.factor(y)) 
        stop("'y' must be a numeric vector")
    if (!((is.vector(y.censored, mode = "numeric") && !is.factor(y.censored)) || 
        is.vector(y.censored, mode = "logical"))) 
        stop("'y.censored' must be a logical or numeric vector")
    if (length(y.censored) != length(y)) 
        stop("'y.censored' must be the same length as 'y'")
    if ((bad.obs.y <- sum(!(ok <- is.finite(y) & is.finite(as.numeric(y.censored))))) > 
        0) {
        y <- y[ok]
        y.censored <- y.censored[ok]
        warning(paste(bad.obs.y, "observations with NA/NaN/Inf in 'y' and/or 'y.censored' removed."))
    }
    if (is.numeric(y.censored)) {
        if (!all(y.censored == 0 | y.censored == 1)) 
            stop(paste("When 'y.censored' is a numeric vector, all values of", 
                "'y.censored' must be 0 (not censored) or 1 (censored)."))
        y.censored <- as.logical(y.censored)
    }
    n.y.cen <- sum(y.censored)
    if (n.x.cen == 0 && n.y.cen == 0) 
        stop(paste("No censored values indicated by 'x.censored'", 
            "and 'y.censored';", "use \n\t\t\tthe function 'wilcox.test' or", 
            "'twoSampleLinearRankTest'"))
    Si.fcn <- function(ni, di, surv.est) {
        k <- length(ni)
        Si <- switch(surv.est, prentice = cumprod((ni - di + 
            1)/(ni + 1)), `kaplan-meier` = cumprod((ni - di)/ni), 
            `peto-peto` = {
                S.hat <- cumprod((ni - di)/ni)
                (S.hat + c(1, S.hat[1:(k - 1)]))/2
            }, altshuler = exp(-cumsum(di/ni)))
    }
    scores.fcn <- function(ni, di, Si, test) {
        k <- length(ni)
        Fi <- 1 - Si
        i <- 1:k
        switch(test, logrank = {
            ci <- -log(Si) - 1
            Ci <- -log(Si)
        }, gehan = {
            ci <- i - ni
            Ci <- i
        }, `tarone-ware` = {
            ci <- i - sqrt(ni)
            Ci <- i
        }, `peto-peto` = {
            ci <- 2 * Fi - 1
            Ci <- Fi
        }, normal.scores.1 = {
            ci <- qnorm(Fi)
            Ci <- dnorm(qnorm(Fi))/Si
        }, normal.scores.2 = {
            ci <- qnorm(Fi)
            Ci <- numeric(k)
            Ci[1] <- -ci[1]/(ni[1] - 1)
            for (j in 2:(k - 1)) Ci[j] <- (ni[j] * Ci[j - 1] - 
                ci[j])/(ni[j] - 1)
            if (ni[k] > 1) Ci[k] <- (ni[k] * Ci[k - 1] - ci[k])/(ni[k] - 
                1)
        }, generalized.sign = {
            ci <- sign(Fi - 0.5)
            Ci <- ifelse(Fi < 0.5, Fi/(1 - Fi), 1)
        })
        list(ci = ci, Ci = Ci)
    }
    cen.levels.x <- sort(unique(x[x.censored]))
    cen.levels.y <- sort(unique(y[y.censored]))
    if (shift.type == "location" && !missing(location.shift.null)) {
        if ((length(location.shift.null) != 1) || !is.finite(location.shift.null)) 
            stop("'location.shift.null' must be a single finite numeric value")
        x <- x - location.shift.null
    }
    if (shift.type == "scale" && !missing(scale.shift.null)) {
        if ((length(scale.shift.null) != 1) || !is.finite(scale.shift.null) || 
            scale.shift.null <= 0) 
            stop("'scale.shift.null' must be a single finite positive numeric value")
        x <- x/scale.shift.null
    }
    max.z <- max(x, y)
    if (censoring.side == "left") {
        x <- (-x) + max.z + 1
        y <- (-y) + max.z + 1
    }
    n.x <- length(x)
    n.y <- length(y)
    N <- n.x + n.y
    z <- c(x, y)
    x.cen <- x[x.censored]
    y.cen <- y[y.censored]
    z.cen <- c(x.cen, y.cen)
    x.no.cen <- x[!x.censored]
    y.no.cen <- y[!y.censored]
    z.no.cen <- c(x.no.cen, y.no.cen)
    n.x.no.cen <- n.x - n.x.cen
    n.y.no.cen <- n.y - n.y.cen
    ti <- sort(z.no.cen)
    k <- length(ti)
    if (variance == "hypergeometric") {
        ti <- unique(ti)
        k <- length(ti)
        ni <- sapply(1:k, function(i, z, ti) sum(ti[i] <= z), 
            z = z, ti = ti)
        n1i <- sapply(1:k, function(i, x, ti) sum(ti[i] <= x), 
            x = x, ti = ti)
        di <- table(z.no.cen)
        d1i <- numeric(k)
        tab <- table(x.no.cen)
        d1i[match(names(tab), as.character(ti))] <- tab
        Si <- Si.fcn(ni = ni, di = di, surv.est = surv.est)
        scores.list <- scores.fcn(ni = ni, di = di, Si = Si, 
            test = test)
        ci <- scores.list$ci
        Ci <- scores.list$Ci
        wi <- ci - Ci
        nu <- sum(wi * (d1i - (di * n1i)/ni))
        r <- n1i/ni
        term <- di * wi^2 * r * (1 - r) * ((ni - di)/(ni - 1))
        if (ni[k] == 1) 
            term[k] <- 0
        var.nu <- sum(term)
    }
    else {
        rle.ti <- rle(ti)
        uncensored.ties <- !all(rle.ti$lengths == 1)
        if (uncensored.ties) {
            diffs <- diff(sort(z))
            delta <- min(diffs[diffs > 0])/(2 * N)
            delta <- min(delta, 1e-08)
            tie.values <- rle.ti$values[rle.ti$lengths > 1]
            tie.values.n <- rle.ti$lengths[rle.ti$lengths > 1]
            n.tie.values <- length(tie.values)
            new.z.no.cen <- z.no.cen
            for (i in 1:n.tie.values) {
                index <- z.no.cen == tie.values[i]
                new.z.no.cen[index] <- z.no.cen[index] - ((tie.values.n[i] - 
                  1):0) * delta
            }
            new.x.no.cen <- new.z.no.cen[1:n.x.no.cen]
            new.y.no.cen <- new.z.no.cen[(n.x.no.cen + 1):k]
            new.ti <- sort(new.z.no.cen)
            new.x <- x
            new.x[!x.censored] <- new.x.no.cen
            new.y <- y
            new.y[!y.censored] <- new.y.no.cen
            new.z <- c(new.x, new.y)
            ni <- sapply(1:k, function(i, z, ti) sum(ti[i] <= 
                z), z = new.z, ti = new.ti)
            n1i <- sapply(1:k, function(i, x, ti) sum(ti[i] <= 
                x), x = new.x, ti = new.ti)
            di <- table(new.z.no.cen)
            d1i <- numeric(k)
            tab <- table(new.x.no.cen)
            d1i[match(names(tab), as.character(new.ti))] <- tab
            Si <- Si.fcn(ni = ni, di = di, surv.est = surv.est)
            scores.list <- scores.fcn(ni = ni, di = di, Si = Si, 
                test = test)
            ci <- scores.list$ci
            Ci <- scores.list$Ci
            ci <- sapply(split(ci, ti), mean)
            Ci <- sapply(split(Ci, ti), mean)
            if (variance == "asymptotic") {
                ai <- cumprod((ni + 1)/(ni + 2))
                ai <- sapply(split(ai, ti), mean)
                Si <- sapply(split(Si, ti), mean)
            }
            ti <- unique(ti)
            k <- length(ti)
            ni <- sapply(1:k, function(i, z, ti) sum(ti[i] <= 
                z), z = z, ti = ti)
            n1i <- sapply(1:k, function(i, x, ti) sum(ti[i] <= 
                x), x = x, ti = ti)
            di <- table(z.no.cen)
            d1i <- numeric(k)
            tab <- table(x.no.cen)
            d1i[match(names(tab), as.character(ti))] <- tab
            ei <- sapply(1:k, function(i, z, ti) sum(ti[i] <= 
                z & z < ti[i + 1]), z = z.cen, ti = c(ti, Inf))
            if (n.x.cen > 0) 
                e1i <- sapply(1:k, function(i, x, ti) sum(ti[i] <= 
                  x & x < ti[i + 1]), x = x.cen, ti = c(ti, Inf))
            else e1i <- rep(0, k)
            nu <- sum(d1i * ci + e1i * Ci)
            var.nu <- switch(variance, permutation = ((n.x * 
                n.y)/(N * (N - 1))) * sum(di * ci^2 + ei * Ci^2), 
                , asymptotic = {
                  bi <- 2 * d1i + e1i
                  dum <- sum(Si * bi) - cumsum(Si * bi)
                  sum(Si * (1 - ai) * bi - (ai - Si) * bi * (Si * 
                    bi + 2 * dum))
                })
        }
        else {
            ni <- sapply(1:k, function(i, z, ti) sum(ti[i] <= 
                z), z = z, ti = ti)
            n1i <- sapply(1:k, function(i, x, ti) sum(ti[i] <= 
                x), x = x, ti = ti)
            di <- table(z.no.cen)
            d1i <- numeric(k)
            tab <- table(x.no.cen)
            d1i[match(names(tab), as.character(ti))] <- tab
            ei <- sapply(1:k, function(i, z, ti) sum(ti[i] <= 
                z & z < ti[i + 1]), z = z.cen, ti = c(ti, Inf))
            if (n.x.cen > 0) 
                e1i <- sapply(1:k, function(i, x, ti) sum(ti[i] <= 
                  x & x < ti[i + 1]), x = x.cen, ti = c(ti, Inf))
            else e1i <- rep(0, k)
            Si <- Si.fcn(ni = ni, di = di, surv.est = surv.est)
            scores.list <- scores.fcn(ni = ni, di = di, Si = Si, 
                test = test)
            ci <- scores.list$ci
            Ci <- scores.list$Ci
            nu <- sum(d1i * ci + e1i * Ci)
            var.nu <- switch(variance, permutation = ((n.x * 
                n.y)/(N * (N - 1))) * sum(di * ci^2 + ei * Ci^2), 
                , asymptotic = {
                  ai <- cumprod((ni + 1)/(ni + 2))
                  bi <- 2 * d1i + e1i
                  dum <- sum(Si * bi) - cumsum(Si * bi)
                  sum(Si * (1 - ai) * bi - (ai - Si) * bi * (Si * 
                    bi + 2 * dum))
                })
        }
    }
    z <- nu/sqrt(var.nu)
    if (censoring.side == "left") {
        nu <- -nu
        z <- -z
    }
    p.value <- switch(alternative, two.sided = 2 * (1 - pnorm(abs(z))), 
        less = pnorm(z), greater = 1 - pnorm(z))
    stat <- c(nu, var.nu, z)
    names(stat) <- c("nu", "var.nu", "z")
    parameters <- NULL
    string <- switch(alternative, two.sided = "!=", less = "<", 
        greater = ">")
    null.value <- "Fx(t)"
    if (shift.type == "location") 
        names(null.value) <- ifelse(missing(location.shift.null), 
            "Fy(t)", paste("Fy(t - ", location.shift.null, ")", 
                sep = ""))
    else {
        names(null.value) <- ifelse(missing(scale.shift.null), 
            "Fy(t)", paste("Fy(t / ", scale.shift.null, ")", 
                sep = ""))
    }
    alternative <- paste(names(null.value), string, "Fx(t) for at least one t")
    sep.string <- paste("\n", space(33), sep = "")
    string0 <- paste(switch(surv.est, prentice = "Prentice", 
        `kaplan-meier` = "Kaplan-Meier", `peto-peto` = "Peto-Peto", 
        altshuler = "Altshuler"), "Survival Estimator")
    string1 <- switch(test, gehan = "Gehan's Test", logrank = "Logrank Test", 
        `tarone-ware` = "Tarone-Ware Test", `peto-peto` = paste("Peto-Peto Test Using", 
            string0, sep = sep.string), normal.scores.1 = paste("Normal Scores Test Using", 
            string0, "Based on Prentice (1978)", sep = sep.string), 
        normal.scores.2 = paste("Normal Scores Test Using", string0, 
            "Based on Prentice and Marek (1979)", sep = sep.string), 
        generalized.sign = "Generalized Sign Test")
    string2 <- switch(variance, permutation = "with Permutation Variance", 
        hypergeometric = "with Hypergeometric Variance", asymptotic = "with Asymptotic Variance")
    method <- paste("Two-Sample Linear Rank Test:", string1, 
        string2, sep = sep.string)
    ret.list <- list(statistic = stat, parameters = parameters, 
        p.value = p.value, estimate = NULL, null.value = null.value, 
        alternative = alternative, method = method, estimation.method = NULL, 
        sample.size = c(nx = n.x, ny = n.y), data.name = data.name, 
        bad.obs = c(x = bad.obs.x, y = bad.obs.y), censoring.side = censoring.side, 
        censoring.name = censoring.name, censoring.levels = list(x = cen.levels.x, 
            y = cen.levels.y), percent.censored = c(x = (100 * 
            n.x.cen)/n.x, y = (100 * n.y.cen)/n.y))
    oldClass(ret.list) <- "htestCensored"
    ret.list
}
