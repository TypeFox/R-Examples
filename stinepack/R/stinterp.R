"stinterp" <-
function (x, y, xout, yp, method = c("scaledstineman", "stineman", 
    "parabola")) 
{
    if (missing(x) || missing(y) || missing(xout)) 
        stop("Wrong number of input arguments, x, y and xout must be specified")
    if (!is.vector(x) || !is.vector(y) || !is.vector(xout) || 
        !is.numeric(x) || !is.numeric(y) || !is.numeric(xout)) 
        stop("x, y and xout must be numeric vectors")
    if (length(x) < 2) 
        stop("x must have 2 or more elements")
    if (length(x) != length(y)) 
        stop("x must have the same number of elements as y")
    if (any(is.na(x)) || any(is.na(x)) || any(is.na(xout))) 
        stop("NAs in x, y or xout are not allowed")
    if (!missing(yp)) {
        if (!is.vector(yp) || !is.numeric(yp)) 
            stop("yp must be a numeric vector")
        if (length(y) != length(yp)) 
            stop("When specified, yp must have the same number of elements as y")
        if (any(is.na(yp)))
           stop("NAs in yp are not allowed")
        if (!missing(method)) 
            stop("Method should not be specified if yp is given")
    }
    dx <- diff(x)
    dy <- diff(y)
    if (any(dx <= 0)) 
        stop("The values of x must strictly increasing")

#calculation of slopes if needed
    if (missing(yp)) {
        yp <- switch(match.arg(method), # this allows for partial argument matching 
               scaledstineman = stinemanSlopes(x,y, scale = TRUE), 
               stineman = stinemanSlopes(x, y, scale = FALSE), 
               parabola = parabolaSlopes(x, y))
    }

# preparations
    m <- length(x)
    m1 <- m - 1
    s <- dy/dx
    k <- length(xout)

    ix <- findInterval(xout, x, rightmost.closed = TRUE)

# For edgepoints allow extrapolation 
# within a tiny range (set by machine precision).

    epx <- 5 * (.Machine$double.eps) * diff(range(x))
    ix[min(x) - epx <= xout & xout <= min(x)] <- 1
    ix[max(x) <= xout & xout <= max(x) + epx] <- m1
    idx <- 1 <= ix & ix <= m1
    ix1 <- ix[idx]
    ix2 <- ix1 + 1

# computation of the interpolant for the three cases dyo1dyo2 ==, > and < 0

    dxo1 <- xout[idx] - x[ix1]
    dxo2 <- xout[idx] - x[ix2]
    y0o <- y[ix1] + s[ix1] * dxo1
    dyo1 <- (yp[ix1] - s[ix1]) * dxo1
    dyo2 <- (yp[ix2] - s[ix1]) * dxo2
    dyo1dyo2 <- dyo1 * dyo2
    yo <- y0o
    if (m > 2 || !missing(yp)) { # linear interpolation is sufficient for m=2 unless slopes are given, then nothing more is done
        id <- dyo1dyo2 > 0
        yo[id] <- y0o[id] + dyo1dyo2[id]/(dyo1[id] + dyo2[id])
        id <- dyo1dyo2 < 0
        yo[id] <- y0o[id] + dyo1dyo2[id] * (dxo1[id] + dxo2[id])/(dyo1[id] - 
            dyo2[id])/((dx[ix1])[id])
    }

# return the results

    yout <- rep(NA, k)
    yout[idx] <- yo
    list(x = xout, y = yout)
}
