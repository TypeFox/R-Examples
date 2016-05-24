
packPotts <- function(x, ncolor) {
    stopifnot(is.matrix(x))
    stopifnot(is.numeric(x))
    ncolor <- as.integer(ncolor)
    storage.mode(x) <- "integer"
    stopifnot(2 <= ncolor && ncolor <= 256)
    stopifnot(all(0 < x & x <= ncolor))
    bitsperpixel <- ceiling(log2(ncolor))
    pixelsperbyte <- floor(8 / bitsperpixel)
    lenraw <- ceiling(length(x) / pixelsperbyte) + 10
    .C("packPotts", x = x, nrow = nrow(x), ncol = ncol(x), ncolor = ncolor,
        lenraw = as.integer(lenraw), raw = raw(lenraw), PACKAGE = "potts")$raw
}

inspectPotts <- function(raw) {
    stopifnot(is.raw(raw))
    foo <- .C("inspectPotts", raw = raw,
        ncolor = integer(1), nrow = integer(1), ncol = integer(1),
        PACKAGE = "potts")
    return(list(ncolor = foo$ncolor, nrow = foo$nrow, ncol = foo$ncol))
}

unpackPotts <- function(raw) {
    stopifnot(is.raw(raw))
    foo <- inspectPotts(raw)
    .C("unpackPotts", raw = raw, lenraw = length(raw),
        ncolor = foo$ncolor, nrow = foo$nrow, ncol = foo$ncol,
        x = matrix(as.integer(0), foo$nrow, foo$ncol), PACKAGE = "potts")$x
}

