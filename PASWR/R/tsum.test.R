tsum.test <-
function(mean.x, s.x = NULL, n.x = NULL, mean.y = NULL, s.y = NULL, n.y = NULL,
    alternative = "two.sided", mu = 0, var.equal = FALSE, conf.level = 0.95)
{
    alt.expanded <- if(!missing(alternative)) char.expand(alternative,
            c("two.sided", "greater", "less"), stop(
            "argument 'alternative' must match one of \"greater\", \"less\", \"two.sided\"."
            )) else alternative
    if(!missing(mu))
        if((length(mu) != 1) || !is.finite(mu))
            stop("argument 'mu' must be a single finite numeric value."
                )
    if(!missing(conf.level))
        if((length(conf.level) != 1) || !is.finite(conf.level) || (
            conf.level <= 0) || (conf.level >= 1))
            stop("argument 'conf.level' must be a single number greater than zero and less than one \n.")
    if(!is.null(mean.x) && is.null(mean.y) && is.null(n.x) && is.null(
        s.x))
        stop("You must enter the value for both s.x and n.x")
    if(is.null(n.x) && !is.null(mean.x) && !is.null(s.x) && is.null(mean.y
        ))
        stop("You must enter the value for n.x")
    if(is.null(s.x) && !is.null(mean.x) && !is.null(n.x) && is.null(mean.y
        ))
        stop("You must enter the value for s.x")
    if(is.null(n.y) && !is.null(mean.x) && !is.null(mean.y) && !is.null(
        s.y) && !is.null(s.x) && !is.null(n.x))
        stop("You must enter the value for n.y")
    if(is.null(n.y) && is.null(n.x) && !is.null(mean.x) && !is.null(mean.y
        ) && !is.null(s.y) && !is.null(s.x))
        stop("You must enter the value for both n.x and n.y")
    if(is.null(s.x) && is.null(s.y) && !is.null(mean.x) && !is.null(mean.y
        ) && !is.null(n.x) && !is.null(n.y))
        stop("You must enter the value for both s.x and s.y")
    if(!is.null(s.x) && is.null(s.y) && !is.null(mean.x) && !is.null(
        mean.y) && !is.null(n.x) && !is.null(n.y))
        stop("You must enter the value for s.y")
    if(is.null(n.y) && is.null(s.y) && !is.null(mean.x) && !is.null(mean.y
        ) && !is.null(s.x) && !is.null(n.x))
        stop("You must enter the value for both s.y and n.y")
    alpha <- 1 - conf.level
    if(is.null(mean.y)) {
        # one-sample t-test.
        # if(var.equal) warning(
        #        "argument 'var.equal' ignored for one-sample test."
        #        )
        conf.int.xbar <- mean.x
        conf.int.s <- sqrt(s.x^2/n.x)
        ret.val <- list(statistic = (conf.int.xbar - mu)/conf.int.s,
            parameters = n.x - 1, estimate = conf.int.xbar,
            null.value = mu, alternative = alt.expanded, method =
            "One-sample t-Test", data.name = c("Summarized x"))
        names(ret.val$estimate) <- "mean of x"
        names(ret.val$null.value) <- "mean"
    }
    else {
        # a two-sample test
        mean.x <- mean.x
        mean.y <- mean.y
        conf.int.xbar <- mean.x - mean.y
        var.x <- s.x^2
        var.y <- s.y^2
        conf.int.s <- if(var.equal) sqrt((((n.x - 1) * var.x + (n.y -
                1) * var.y) * (1/n.x + 1/n.y))/(n.x + n.y -
                2)) else sqrt((var.x/n.x) + (var.y/n.y))
        ret.val <- c(if(var.equal) list(method =
                "Standard Two-Sample t-Test", parameters = n.x +
                n.y - 2) else list(method =
                "Welch Modified Two-Sample t-Test", parameters
                 = {
                const <- 1/(1 + (n.x * var.y)/(n.y * var.x))
                1/((const^2)/(n.x - 1) + ((1 - const)^2)/
                    (n.y - 1))
            }
            ), list(statistic = (conf.int.xbar - mu)/conf.int.s,
            estimate = c(mean.x, mean.y), null.value = mu,
            alternative = alt.expanded,
            data.name = paste("Summarized ", deparse(substitute (x)), " and ", deparse(substitute(y)), sep = "")))
        names(ret.val$estimate) <- c("mean of x", "mean of y")
        names(ret.val$null.value) <- "difference in means"
    }
    ret.val <- c(ret.val, switch(alt.expanded,
        two.sided = {
            conf.int.hw <- qt((1 - alpha/2), ret.val$parameters) *
                conf.int.s
            list(p.value = 2 * pt( - abs(ret.val$statistic),
                ret.val$parameters), conf.int = c(
                conf.int.xbar - conf.int.hw, conf.int.xbar +
                conf.int.hw))
        }
        ,
        greater = {
            list(p.value = 1 - pt(ret.val$statistic, ret.val$
                parameters), conf.int = c(conf.int.xbar - qt(
                (1 - alpha), ret.val$parameters) * conf.int.s,
                Inf))
        }
        ,
        less = {
            list(p.value = pt(ret.val$statistic, ret.val$
                parameters), conf.int = c(-Inf, conf.int.xbar +
                qt((1 - alpha), ret.val$parameters) *
                conf.int.s))
        }
        ))
    names(ret.val$statistic) <- "t"
    names(ret.val$parameters) <- "df"
    attr(ret.val$conf.int, "conf.level") <- conf.level
    ret.val <- ret.val[c("statistic", "parameters", "p.value", "conf.int",
        "estimate", "null.value", "alternative", "method", "data.name"
        )]
    oldClass(ret.val) <- "htest"
    return(ret.val)
}

