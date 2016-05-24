if(! exists(".bincode", envir = .BaseNamespaceEnv))
    .bincode <- function(v, breaks, ...) {
        .C("bincode", as.double(v), length(v), as.double(breaks),
           length(breaks), code = integer(length(v)), as.logical(TRUE),
           as.logical(TRUE), nok = TRUE, NAOK = TRUE, DUP = FALSE,
           PACKAGE = "base")$code
    }

image3d <- function (v, x = 1:dim(v)[1], y = 1:dim(v)[2], z = 1:dim(v)[3],
                     vlim = quantile(v, c(.9, 1),na.rm=TRUE),
                     col = heat.colors(256),
                     alpha.power = 2,
                     alpha = ((1:length(col))/ length(col))^alpha.power,
                     breaks,sprites = TRUE, jitter = FALSE,
                     radius = min(diff(x), diff(y), diff(z)),
                     add = FALSE,...)
{
    loadRGL()
    if (!is.array(v) && length(dim(v)) != 3)
        stop("'v' must be a 3D array")
    nx <- dim(v)[1]
    ny <- dim(v)[2]
    nz <- dim(v)[3]
    if (length(x) != nx || length(y) != ny || length(z) != nz)
        stop("dimensions of v do not match x, y, or z")
    if (missing(breaks)) {
        nc <- length(col)
        if (any(!is.finite(vlim)) || diff(vlim) < 0)
            stop("invalid v limits")
        if (diff(vlim) == 0)
            vlim <- if (vlim[1] == 0)
                c(-1, 1)
            else vlim[1] + c(-0.4, 0.4) * abs(vlim[1])
        v <- (v - vlim[1])/diff(vlim)
        vi <- floor((nc - 1e-05) * v + 1e-07)
        vi[vi < 0 | vi >= nc] <- NA
        if (length(alpha) == 1)
            alpha <- rep(alpha, nc)
        else if (length(alpha) != nc)
            stop("number of colors and alpha levels must be identical")
    }
    else {
        if (length(breaks) != length(col) + 1)
            stop("must have one more break than colour")
        if (length(breaks) != length(alpha) + 1)
            stop("must have one more break than alpha levels")
        if (any(!is.finite(breaks)))
            stop("breaks must all be finite")
        vi <- .bincode(v, breaks, TRUE, TRUE) - 1
    }
    if (!add)
        clear3d()
    i <- which(is.finite(vi))
    xi <- x[as.integer((i - 1) %% nx + 1)]
    yi <- y[as.integer(((i - 1) %/% nx) %% ny + 1)]
    zi <- z[as.integer((i - 1) %/% (nx * ny) + 1)]
    vi <- vi[i] + 1
    if (jitter) {
      ni <- length(i)
      xi <- xi + runif(ni, max = min(diff(x)))
      yi <- yi + runif(ni, max = min(diff(y)))
      zi <- zi + runif(ni, max = min(diff(z)))
    }
    if (sprites) {
        texture <- system.file("textures/particle.png", package="rgl")
        sprites3d(xi, yi, zi, color = col[vi], alpha = alpha[vi],
                    lit=FALSE, radius = radius, textype="alpha",
                    texture = texture, ...)
    }
    else points3d(xi, yi, zi, color = col[vi], alpha = alpha[vi], ...)
}
