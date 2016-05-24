zsum.test <-
function(mean.x, sigma.x = NULL, n.x = NULL, mean.y = NULL, sigma.y = NULL,
    n.y = NULL, alternative = "two.sided", mu = 0, conf.level = 0.95)
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
    if(!is.null(mean.x) && is.null(mean.y) && is.null(n.x) && is.null(
        sigma.x))
        stop("You must enter the value for both sigma.x and n.x")
    if(is.null(n.x) && !is.null(mean.x) && !is.null(sigma.x) && is.null(
        mean.y))
        stop("You must enter the value for n.x")
    if(is.null(sigma.x) && !is.null(mean.x) && !is.null(n.x) && is.null(
        mean.y))
        stop("You must enter the value for sigma.x")
    if(is.null(n.y) && !is.null(mean.x) && !is.null(mean.y) && !is.null(
        sigma.y) && !is.null(sigma.x) && !is.null(n.x))
        stop("You must enter the value for n.y")
    if(is.null(n.y) && is.null(n.x) && !is.null(mean.x) && !is.null(mean.y
        ) && !is.null(sigma.y) && !is.null(sigma.x))
        stop("You must enter the value for both n.x and n.y")
    if(is.null(sigma.x) && is.null(sigma.y) && !is.null(mean.x) && !
        is.null(mean.y) && !is.null(n.x) && !is.null(n.y))
        stop("You must enter the value for both sigma.x and sigma.y")
    if(!is.null(sigma.x) && is.null(sigma.y) && !is.null(mean.x) && !
        is.null(mean.y) && !is.null(n.x) && !is.null(n.y))
        stop("You must enter the value for sigma.y")
    if(is.null(n.y) && is.null(sigma.y) && !is.null(mean.x) && !is.null(
        mean.y) && !is.null(sigma.x) && !is.null(n.x))
        stop("You must enter the value for both sigma.y and n.y")
    if(!missing(conf.level))
        if(length(conf.level) != 1 || is.na(conf.level) || conf.level <
            0 || conf.level > 1)
            stop("conf.level must be a number between 0 and 1")
    if(!is.null(mean.y)) {
        dname <- c("Summarized x and y")
    }
    else {
        dname <- c("Summarized x")
    }
    n.x
    if(n.x <= 1)
        stop("not enough x observations")
    estimate <- mean.x
    if(is.null(mean.y)) {
        stderr <- sigma.x/sqrt(n.x)
        zobs <- (mean.x - mu)/stderr
        method <- c("One-sample z-Test")
        names(estimate) <- c("mean of x")
    }
    else {
        n.y
        if(n.y <= 1)
            stop("not enough y observations")
        method <- c("Two-sample z-Test")
        estimate <- c(mean.x, mean.y)
        names(estimate) <- c("mean of x", "mean of y")
        stderr <- sqrt(((sigma.x^2)/n.x) + ((sigma.y^2)/n.y))
        zobs <- (mean.x - mean.y - mu)/stderr
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
    if(!is.null(mean.y))
        names(mu) <- "difference in means"
    else names(mu) <- "mean"
    attr(cint, "conf.level") <- conf.level
    rval <- list(statistic = zobs, p.value = pval, conf.int = cint,
        estimate = estimate, null.value = mu, alternative =
        alternative, method = method, data.name = dname)
    attr(rval, "class") <- "htest"
    return(rval)
}

