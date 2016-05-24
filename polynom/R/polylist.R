polylist <-
function(...)
    .polylist_from_list(list(...))

.polylist_from_list <-
function(x)
    structure(lapply(x, as.polynomial), class = "polylist")

is.polylist <-
function(x)
    inherits(x, "polylist")

as.polylist <-
function(x)
{
    if(is.polylist(x)) x
    else if(is.list(x)) .polylist_from_list(x)
    else polylist(x)
}

deriv.polylist <-
function(expr, ...) 
    structure(lapply(expr, deriv), class = class(expr))

integral.polylist <-
function(expr, ...)
{
    result <- lapply(expr, integral, ...)
    if (length(result) > 0 && is.polynomial(result[[1]]))
        class(result) <- class(expr)
    result
}

plot.polylist <-
function(x, xlim = 0:1, ylim = range(Px), type = "l", len = 100, ...)
{
    p <- x                              # generic/method
    if(missing(xlim)) {
        ## try to cover the "interesting" region
        xlim <- range(Re(unlist(lapply(p, summary.polynomial))))
    }
    if(any(is.na(xlim))) {
        warning("summary of polynomial fails. Using nominal xlim")
        xlim <- 0:1
    }
    if(diff(xlim) == 0)
        xlim <- xlim + c(-1, 1)/2
    if(length(xlim) > 2)
        x <- xlim
    else {
        eps <- diff(xlim)/100
        xlim <- xlim + c( - eps, eps)
        x <- seq(xlim[1], xlim[2], len = len)
    }
    Px <- unlist(lapply(p, predict.polynomial, x))
    if(!missing(ylim))
        Px[Px < ylim[1]] <- Px[Px > ylim[2]] <- NA
    plot(cbind(x, Px), xlab = "x", ylab = "P(x)", type = "n",
         xlim = xlim, ylim = ylim, ...)
    for(i in seq(along = p))
        lines(p[[i]], lty = i)
    invisible()
}

print.polylist <-
function(x, ...)
{
    cat("List of polynomials:\n")
    y <- x
    x <- unclass(x)
    NextMethod()
    invisible(y)
}

c.polylist <-
function(..., recursive = FALSE)
    .polylist_from_list(unlist(lapply(list(...), as.polylist),
                               recursive = FALSE))

"[.polylist" <-
function(x, i)
    .polylist_from_list(NextMethod("["))

rep.polylist <-
function(x, times, ...)
    .polylist_from_list(NextMethod("rep"))

unique.polylist <-
function(x, incomparables = FALSE, ...)
    .polylist_from_list(NextMethod("unique"))

Summary.polylist <-
function(..., na.rm = FALSE)
{
    ok <- switch(.Generic, sum = , prod = TRUE, FALSE)
    if(!ok)
        stop(gettextf("Generic '%s' not defined for \"%s\" objects.",
                      .Generic, .Class))
    switch(.Generic,
           "sum" = accumulate("+", as.polynomial(0), c(...)),
           "prod" = accumulate("*", as.polynomial(1), c(...)))
}

