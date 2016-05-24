##' @title Replicate and cbind all block-level permutations
##' @param x a list whose compontents are the set of all permutations
##' at the block level
##' @return a matrix
##' @author Gavin L. Simpson
`cbindAllPerms` <- function(x) {
    nb <- length(x) ## number of blocks

    ## allPerms has first block varying slowest, but expand.grid has
    ## first factor varying fastest. Hence we reverse `x` here, and
    ## also reverse `out` back again later
    x <- rev(x)

    ## prepares nb seqence vectors 1:`obs in block` for expand.grid
    rowind <- do.call(expand.grid,
                      lapply(x, function(i) seq_len(nrow(i))))

    ## index elements of x using the row indices - gives a list to cbind
    ## next. sapply() over-simplifies to wrong dimensions so not used.
    ## Note: the lapply() result is reversed with `rev` to undo the reversing
    ## of `x` earlier; ensures everything is correct block order.
    out <- rev(lapply(seq_len(nb), function(i, m, ind) m[[i]][ind[, i] ,],
                      m = x, ind = rowind))
    do.call(cbind, out) ## return
}
