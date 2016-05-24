setClass("smoothbin",
         representation("hexbin", wts="numeric"))

setGeneric("hsmooth", function(bin, wts) standardGeneric("hsmooth"))

smooth.hexbin <- function(bin, wts = c(48, 4, 1))
{
    if(!is(bin,"hexbin"))
        stop("first argument must be a hexbin object")
    cell <- bin@cell - 1
    n <- as.integer(length(cell))
    cnt <- bin@count
    xbins <- bin@xbins
    bdim <- bin@dimen
    row <- bdim[1]
    col <- bdim[2]
    ncol <- col + 4
    nrow <- row + 4
    nmax <- ncol * nrow
    sm <- rep.int(0:0, nmax)
    nr <- cell %/% col + 2
    nc <- cell %% col + 3
    pad <- rep.int(0:0, nmax - n)
    cell <- c(nr * ncol + nc, pad)
    cnt <- c(cnt, pad)
    ans <- .Fortran("hsm",
                    cell = as.integer(cell),
                    cnt  = as.integer(cnt),
                    n    = n,
                    nmax = as.integer(nmax),
                    sm   = as.integer(sm),
                    ncol = as.integer(ncol),
                    wts  = as.integer(wts),
                    PACKAGE = "hexbin")[c("cell","cnt","n")]
    n <- ans$n
    length(ans$cell) <- length(ans$cnt) <- n
    bin@xbins <- xbins + 4
    bin@xbnds <- addBit(bin@xbnds, f = 2/xbins)
    bin@ybnds <- addBit(bin@ybnds, f = 2/xbins)
    bin@dimen <- c(nrow, ncol)
    bin@cell <- ans$cell
    bin@count <- ans$cnt
    new("smoothbin", bin, wts=wts)
}

setMethod("hsmooth", "hexbin", smooth.hexbin)
