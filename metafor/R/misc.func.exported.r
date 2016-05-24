`[.escalc` <-
function (x, ...) 
{
    dat <- NextMethod("[")
    yi.names <- attr(x, "yi.names")
    if (!is.null(yi.names)) {
        for (i in 1:length(yi.names)) {
            if (!is.element(yi.names[i], names(dat))) {
                next
            }
            else {
                eval(parse(text = paste0("attr(dat$", yi.names[i], 
                  ", 'measure') <- attr(x$", yi.names[i], ", 'measure')")))
            }
        }
    }
    attr(dat, "yi.names") <- attr(x, "yi.names")
    attr(dat, "vi.names") <- attr(x, "vi.names")
    attr(dat, "sei.names") <- attr(x, "sei.names")
    attr(dat, "zi.names") <- attr(x, "zi.names")
    attr(dat, "ci.lb.names") <- attr(x, "ci.lb.names")
    attr(dat, "ci.ub.names") <- attr(x, "ci.ub.names")
    if (!is.null(attr(x, "digits"))) 
        attr(dat, "digits") <- attr(x, "digits")
    return(dat)
}
`[.list.rma` <-
function (x, i, ...) 
{
    out <- x
    attr(out, "class") <- NULL
    slab.pos <- which(names(out) == "slab")
    if (!missing(i)) 
        out[seq_len(slab.pos - 1)] <- lapply(out[seq_len(slab.pos - 
            1)], function(r) if (class(r) == "matrix") 
            r[i, ]
        else r[i])
    if (length(out[[1]]) == 0L) 
        return(NULL)
    out$slab <- x$slab[i]
    if (anyNA(out$slab)) 
        return(NULL)
    out$digits <- x$digits
    out$transf <- x$transf
    out$method <- x$method
    class(out) <- "list.rma"
    return(out)
}
cbind.escalc <-
function (..., deparse.level = 1) 
{
    dat <- data.frame(..., check.names = FALSE)
    arguments <- list(...)
    yi.names <- NULL
    vi.names <- NULL
    sei.names <- NULL
    zi.names <- NULL
    ci.lb.names <- NULL
    ci.ub.names <- NULL
    digits <- NULL
    for (arg in arguments) {
        yi.names <- c(attr(arg, "yi.names"), yi.names)
        vi.names <- c(attr(arg, "vi.names"), vi.names)
        sei.names <- c(attr(arg, "sei.names"), sei.names)
        zi.names <- c(attr(arg, "zi.names"), zi.names)
        ci.lb.names <- c(attr(arg, "ci.lb.names"), ci.lb.names)
        ci.ub.names <- c(attr(arg, "ci.ub.names"), ci.ub.names)
        digits <- c(attr(arg, "digits"), digits)
    }
    attr(dat, "yi.names") <- unique(yi.names)
    attr(dat, "vi.names") <- unique(vi.names)
    attr(dat, "sei.names") <- unique(sei.names)
    attr(dat, "zi.names") <- unique(zi.names)
    attr(dat, "ci.lb.names") <- unique(ci.lb.names)
    attr(dat, "ci.ub.names") <- unique(ci.ub.names)
    attr(dat, "digits") <- digits[1]
    class(dat) <- c("escalc", "data.frame")
    return(dat)
}
replmiss <-
function (x, y) 
{
    if (length(x) == 0) 
        return(x)
    if (length(y) == 1) 
        y <- rep(y, length(x))
    if (length(x) != length(y)) 
        stop("Length of 'x' and 'y' do not match.")
    x[is.na(x)] <- y[is.na(x)]
    return(x)
}
