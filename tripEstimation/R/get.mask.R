"get.mask" <-
function (masks, k)
{
    list(x = masks$x, y = masks$y, z = bits(masks$z[, , (k%/%31) + 1], (k - 1)%%31))
}

"mkMaskObject" <- function(xs, ys, nsegs) {
    nblocks <- (nsegs%/%31) + 1
    nlon <- length(xs)
    nlat <- length(ys)
    list(x = xs, y = ys, z = array(integer(nblocks * nlat * nlon), c(nlon, nlat, nblocks)))
}


"set.mask<-" <- function(object, segment, value) {
    k <- segment
    blk <- (k%/%31) + 1
    bit <- (k - 1)%%31
    bits(object$z[,,blk], bit) <- value
    object
}
