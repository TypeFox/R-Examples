`checkModifications` <-
function (fun, Length = 100) 
{
    p <- seq(.Machine$double.neg.eps, 1 - 1e-10, length = Length)
    te <- fun(p)
    if (!is.list(te)) 
        stop("The result should be a list of length two.")
    if (length(te) != 2) 
        stop("The result should be a list of length two.")
    if (any(is.na(match(names(te), c("ar", "at"))))) 
        stop("The result should be a list with elements 'ar' and'at'.")
    if (length(te$ar) != Length) 
        stop("'ar' should be of the same length as 'p'")
    if (length(te$at) != Length) 
        stop("'at' should be of the same length as 'p'")
    if (any(te$ar >= te$at)) 
        stop("'ar' cannot take larger values than 'at'")
    if (any(te$ar < 0)) 
        stop("'ar' cannot be negative")
    if (any(te$ar < 0)) 
        stop("'at' cannot be negative")
    plot(p, te$at, ylim = c(0, 10), type = "l")
    points(p, te$ar, type = "l", col = "grey")
    drop(TRUE)
}
`modifications` <-
function (family, pl = FALSE) 
{
    if (is.character(family)) 
        family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family)) 
        family <- family()
    distr.link <- paste(family$family, family$link[1], sep = ".")
    distr.link <- gsub(pattern = "[(]", x = distr.link, replacement = ".")
    distr.link <- gsub(pattern = "[)]", x = distr.link, replacement = ".")
    if (pl) {
        out <- switch(distr.link, binomial.logit = function(p) {
            etas <- family$linkfun(p)
            list(ar = 0.5 * p/p, at = 1 * p/p)
        }, binomial.probit = function(p) {
            etas <- family$linkfun(p)
            list(ar = p * (1 - etas * (etas < 0)/dnorm(etas)), 
                at = etas * ((etas >= 0) * (1 - p) - (etas < 
                  0) * p)/dnorm(etas) + 0.5/p)
        }, binomial.cloglog = function(p) {
            etas <- family$linkfun(p)
            list(ar = -p/log(1 - p), at = 0.5/p)
        }, binomial.cauchit = function(p) {
            etas <- family$linkfun(p)
            list(ar = -2 * pi * etas * p * (etas < 0) + (p - 
                0.5) * (etas >= 0) + p, at = 2 * pi * etas * 
                ((etas >= 0) - p) - (p - 0.5)/p * (etas < 0) + 
                1)
        }, NULL)
        if (is.null(out)) 
            out <- match.fun("mpl.custom.family")
    }
    else {
        out <- switch(distr.link, binomial.logit = function(p) {
            etas <- family$linkfun(p)
            list(ar = 0.5 * p/p, at = 1 * p/p)
        }, binomial.probit = function(p) {
            etas <- family$linkfun(p)
            list(ar = -0.5 * p * etas * (etas < 0)/dnorm(etas) + 
                p, at = 0.5 * etas * ((etas >= 0) - p)/dnorm(etas) + 
                1)
        }, binomial.cloglog = function(p) {
            etas <- family$linkfun(p)
            list(ar = -0.5 * p/log(1 - p) + p, at = 0.5 * p/p + 
                1)
        }, binomial.cauchit = function(p) {
            etas <- family$linkfun(p)
            list(ar = -pi * etas * p * (etas < 0) + p, at = pi * 
                etas * ((etas >= 0) - p) + 1)
        }, NULL)
        if (is.null(out)) 
            out <- match.fun("br.custom.family")
    }
    out
}
