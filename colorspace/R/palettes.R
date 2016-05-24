rainbow_hcl <- function(n, c = 50, l = 70, start = 0, end = 360 * (n - 1)/n,
                        gamma = NULL, fixup = TRUE, alpha = 1, ...)
{
    if (!is.null(gamma))
        warning("'gamma' is deprecated and has no effect")
    if(n < 1L) return(character(0L))
    rval <- hex(polarLUV(L = l, C = c, H = seq(start, end, length = n)),
                fixup = fixup, ...)

    if(!missing(alpha)) {
        alpha <- pmax(pmin(alpha, 1), 0)
	alpha <- format(as.hexmode(round(alpha * 255 + 0.0001)),
	                width = 2L, upper.case = TRUE)
        rval <- paste(rval, alpha, sep = "")
    }

    return(rval)
}

diverge_hcl <- function(n, h = c(260, 0), c = 80, l = c(30, 90), power = 1.5,
                        gamma = NULL, fixup = TRUE, alpha = 1, ...)
{
    if (!is.null(gamma))
        warning("'gamma' is deprecated and has no effect")
    if(n < 1L) return(character(0L))
    h <- rep(h, length.out = 2L)
    c <- c[1L]
    l <- rep(l, length.out = 2L)
    power <- rep(power, length.out = 2L)
    rval <- seq(1, -1, length = n)
    rval <- hex(polarLUV(L = l[2L] - diff(l) * abs(rval)^power[2L],
                         C = c * abs(rval)^power[1L],
                         H = ifelse(rval > 0, h[1L], h[2L])),
                fixup = fixup, ...)

    if(!missing(alpha)) {
        alpha <- pmax(pmin(alpha, 1), 0)
	alpha <- format(as.hexmode(round(alpha * 255 + 0.0001)),
	                width = 2L, upper.case = TRUE)
        rval <- paste(rval, alpha, sep = "")
    }

    return(rval)
}

diverge_hsv <- function(n, h = c(240, 0), s = 1, v = 1, power = 1,
                        gamma = NULL, fixup = TRUE, alpha = 1, ...)
{
    if (!is.null(gamma))
        warning("'gamma' is deprecated and has no effect")
    if(n < 1L) return(character(0L))
    h <- rep(h, length.out = 2L)
    s <- s[1L]
    v <- v[1L]
    power <- power[1L]
    rval <- seq(-s, s, length = n)
    rval <- hex(as(HSV(H = ifelse(rval > 0, h[2L], h[1L]),
                       S = abs(rval)^power, V = v, ...), "RGB"),
                fixup = fixup, ...)

    if(!missing(alpha)) {
        alpha <- pmax(pmin(alpha, 1), 0)
	alpha <- format(as.hexmode(round(alpha * 255 + 0.0001)),
	                width = 2L, upper.case = TRUE)
        rval <- paste(rval, alpha, sep = "")
    }

    return(rval)
}

sequential_hcl <- function(n, h = 260, c. = c(80, 0), l = c(30, 90),
                           power = 1.5, gamma = NULL, fixup = TRUE, alpha = 1, ...)
{
    if (!is.null(gamma))
        warning("'gamma' is deprecated and has no effect")
    if(n < 1L) return(character(0L))
    c <- rep(c., length.out = 2L)
    l <- rep(l, length.out = 2L)
    power <- rep(power, length.out = 2L)
    rval <- seq(1, 0, length = n)
    rval <- hex(polarLUV(L = l[2L] - diff(l) * rval^power[2L],
                         C = c[2L] - diff(c) * rval^power[1L],
                         H = h[1L]),
                fixup = fixup, ...)

    if(!missing(alpha)) {
        alpha <- pmax(pmin(alpha, 1), 0)
	alpha <- format(as.hexmode(round(alpha * 255 + 0.0001)),
	                width = 2L, upper.case = TRUE)
        rval <- paste(rval, alpha, sep = "")
    }

    return(rval)
}

heat_hcl <- function(n, h = c(0, 90), c. = c(100, 30), l = c(50, 90),
                     power = c(1/5, 1), gamma = NULL, fixup = TRUE, alpha = 1, ...)
{
    if (!is.null(gamma))
        warning("'gamma' is deprecated and has no effect")
    if(n < 1L) return(character(0L))
    h <- rep(h, length.out = 2L)
    c <- rep(c., length.out = 2L)
    l <- rep(l, length.out = 2L)
    power <- rep(power, length.out = 2L)
    rval <- seq(1, 0, length = n)
    rval <- hex(polarLUV(L = l[2L] - diff(l) * rval^power[2L],
                         C = c[2L] - diff(c) * rval^power[1L],
                         H = h[2L] - diff(h) * rval),
                fixup = fixup, ...)

    if(!missing(alpha)) {
        alpha <- pmax(pmin(alpha, 1), 0)
	alpha <- format(as.hexmode(round(alpha * 255 + 0.0001)),
	                width = 2L, upper.case = TRUE)
        rval <- paste(rval, alpha, sep = "")
    }

    return(rval)
}

terrain_hcl <- function(n, h = c(130, 0), c. = c(80, 0), l = c(60, 95),
                        power = c(1/10, 1), gamma = NULL, fixup = TRUE, alpha = 1, ...)
{
    if (!is.null(gamma))
        warning("'gamma' is deprecated and has no effect")
    rval <- heat_hcl(n, h = h, c. = c., l = l, power = power,
                     fixup = fixup, ...)

    if(!missing(alpha)) {
        alpha <- pmax(pmin(alpha, 1), 0)
	alpha <- format(as.hexmode(round(alpha * 255 + 0.0001)),
	                width = 2L, upper.case = TRUE)
        rval <- paste(rval, alpha, sep = "")
    }

    return(rval)
}
