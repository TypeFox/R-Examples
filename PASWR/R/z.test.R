z.test <-
function(x, y = NULL, alternative = "two.sided", mu = 0, sigma.x = NULL,
    sigma.y = NULL, conf.level = 0.95)
{
    choices <- c("two.sided", "greater", "less")
    alt <- pmatch(alternative, choices)
    alternative <- choices[alt]
    if(length(alternative) > 1 || is.na(alternative))
        stop("alternative must be one \"greater\", \"less\", \"two.sided\""
            )
    if(!missing(mu))
        if(length(mu) != 1 || is.na(mu))
            stop("mu must be a single number")
    if(is.null(sigma.x) && !is.null(x) && is.null(y))
        stop("You must enter the value for sigma.x")
    if(!is.null(y) && is.null(sigma.y) || is.null(sigma.x))
        stop("You must enter values for both sigma.x and sigma.y")
    if(!missing(conf.level))
        if(length(conf.level) != 1 || is.na(conf.level) || conf.level <
            0 || conf.level > 1)
            stop("conf.level must be a number between 0 and 1")
    if(!is.null(y)) {
        dname <- paste(deparse(substitute(x)), "and", paste(deparse(
            substitute(y))))
    }
    else {
        dname <- deparse(substitute(x))
    }
    nx <- length(x)
    if(nx <= 2)
        stop("not enough x observations")
    mx <- mean(x)
    estimate <- mx
    if(is.null(y)) {
        stderr <- sigma.x/sqrt(nx)
        zobs <- (mx - mu)/stderr
        method <- c("One-sample z-Test")
        names(estimate) <- c("mean of x")
    }
    else {
        ny <- length(y)
        if(ny <= 2)
            stop("not enough y observations")
        my <- mean(y)
        method <- c("Two-sample z-Test")
        estimate <- c(mx, my)
        names(estimate) <- c("mean of x", "mean of y")
        stderr <- sqrt(((sigma.x^2)/nx) + ((sigma.y^2)/ny))
        zobs <- (mx - my - mu)/stderr
    }
    if(alternative == "less") {
        pval <- pnorm(zobs)
        cint <- c(-Inf, zobs * stderr + qnorm(conf.level) * stderr)
    }
    else if(alternative == "greater") {
        pval <- 1 - pnorm(zobs)
        cint <- c(zobs * stderr - qnorm(conf.level) * stderr, Inf)
    }
    else {
        pval <- 2 * pnorm( - abs(zobs))
        alpha <- 1 - conf.level
        cint <- c(zobs * stderr - qnorm((1 - alpha/2)) * stderr, zobs *
            stderr + qnorm((1 - alpha/2)) * stderr)
    }
    cint <- cint + mu
    names(zobs) <- "z"
    if(!is.null(y))
        names(mu) <- "difference in means"
    else names(mu) <- "mean"
    attr(cint, "conf.level") <- conf.level
    rval <- list(statistic = zobs, p.value = pval, conf.int = cint,
        estimate = estimate, null.value = mu, alternative =
        alternative, method = method, data.name = dname)
    attr(rval, "class") <- "htest"
    return(rval)
}

