`allGrid` <- function(n, nperms, nr, nc, mirror, constant)
{
    v <- seq_len(n)
    X <- matrix(nrow = nperms, ncol = n)
    idx <- 1
    ## ncol == 2 is special case
    if(nc == 2) {
        X <- allSeries(n, nperms = nperms, mirror = mirror)
    } else {
        for(i in seq_len(nr)) {
            for(j in seq_len(nc)) {
                ir <- seq(i, length = nr)%%nr
                ic <- seq(j, length = nc)%%nc
                ## block 1 - no reversals
                X[idx, ] <- rep(ic, each = nr) * nr +
                    rep(ir, len = nr * nc) + 1
                if(mirror) {
                    ## block 2 - rev rows but not columns
                    X[idx + n, ] <- rep(ic, each = nr) * nr +
                        rep(rev(ir), len = nr * nc) + 1
                    ## block 3 - rev columns but not rows
                    X[idx + (2*n), ] <- rep(rev(ic), each = nr) *
                        nr + rep(ir, len = nr * nc) + 1
                }
                idx <- idx + 1
            }
        }
        if(mirror) {
            ## rev columns and rows
            ## no calculations, just rev cols of block 1
            v <- seq_len(n)
            X[((3*n)+1):(4*n), ] <- X[v, rev(v)]
        }
    }
    X
}
