gevmle <- function(x, ..., method = "Nelder"){

    n <- length(x)
    param <- c("loc", "scale", "shape")

    nlgev <- function(loc, scale, shape)
        -.C("gevlik", as.double(x), as.integer(n),
            as.double(loc), as.double(scale),
            as.double(shape), dns = double(1),
            PACKAGE = "SpatialExtremes", NAOK = TRUE)$dns

    start <- c(loc = 0, scale = 0, shape = 0)
    start["scale"] <- sqrt(6 * var(x, na.rm = TRUE)) / pi
    start["loc"] <- mean(x, na.rm = TRUE) - 0.58 * start["scale"]

    start <- start[!(param %in% names(list(...)))]

    nm <- names(start)
    l <- length(nm)
    f <- formals(nlgev)
    names(f) <- param
    m <- match(nm, param)

    formals(nlgev) <- c(f[m], f[-m])
    nllh <- function(p, ...) nlgev(p, ...)

    if(l > 1)
        body(nllh) <- parse(text = paste("nlgev(", paste("p[",1:l,
                            "]", collapse = ", "), ", ...)"))

    fixed.param <- list(...)[names(list(...)) %in% param]

    if(any(!(param %in% c(nm,names(fixed.param)))))
        stop("unspecified parameters")

    opt <- optim(start, nllh, ..., method = method)
    param <- c(opt$par, unlist(fixed.param))
    names(param) <- c("loc", "scale", "shape")

    return(param)
}


gpdmle <- function(x, threshold, ..., method = "Nelder"){

    param <- c("scale", "shape")

    nlgpd <- function(scale, shape)
        -.C("gpdlik", as.double(exceed), as.integer(nat), as.double(threshold),
            as.double(scale), as.double(shape), dns = double(1),
            PACKAGE = "SpatialExtremes", NAOK = TRUE)$dns

    high <- (x > threshold) & !is.na(x)
    exceed <- as.double(x[high])
    nat <- as.integer(length(exceed))

    start <- c(scale = mean(exceed, na.rm = TRUE) - min(threshold), shape = 0.0)

    start <- start[!(param %in% names(list(...)))]

    nm <- names(start)
    l <- length(nm)
    f <- formals(nlgpd)
    names(f) <- param
    m <- match(nm, param)

    formals(nlgpd) <- c(f[m], f[-m])
    nllh <- function(p, ...) nlgpd(p, ...)

    if(l > 1)
        body(nllh) <- parse(text = paste("nlgpd(", paste("p[",1:l,
                            "]", collapse = ", "), ", ...)"))

    fixed.param <- list(...)[names(list(...)) %in% param]

    if(any(!(param %in% c(nm,names(fixed.param)))))
        stop("unspecified parameters")

    opt <- optim(start, nllh, ..., method = method)
    param <- c(opt$par, unlist(fixed.param))
    names(param) <- c("scale", "shape")

    return(param)
}
