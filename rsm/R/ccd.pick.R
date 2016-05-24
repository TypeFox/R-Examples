### Find good ccds

ccd.pick = function(k, n.c=2^k, n0.c=1:10, blks.c=1, n0.s=1:10, 
    bbr.c=1, wbr.s=1, bbr.s=1, best=10, 
    sortby=c("agreement","N"), restrict)
{
    grid = expand.grid (n.c=n.c, n0.c=n0.c, blks.c=blks.c, n.s=NA, n0.s=n0.s, bbr.c=bbr.c, wbr.s=wbr.s, bbr.s=bbr.s)
    grid$n.s = 2 * k * grid$wbr.s
    grid$N = with(grid, blks.c * bbr.c * (n.c + n0.c) + bbr.s * (2 * k * wbr.s + n0.s))
    grid$alpha.rot = (with(grid, n.c * blks.c * bbr.c / (wbr.s * bbr.s))) ^ .25
    num = with(grid, n.c * (2 * k * wbr.s + n0.s))
    den = with(grid, 2 * wbr.s * (n.c + n0.c))
    grid$alpha.orth = sqrt(num / den)
    agreement = with(grid, abs(log(alpha.rot / alpha.orth)))
    
    # remove combinations that don't have enough d.f.
    extra.df = with(grid, n.c * blks.c - (k * (k + 1)/2 + blks.c))
    grid = grid[extra.df >= 0, ]
    
    if (!missing(restrict))
        for (restr in restrict) {
            r = with(grid, eval(parse(text=restr)))
            grid = grid[r,]
        }
    
    # regenerate 'agreement' in case # rows changed
    agreement = with(grid, abs(log(alpha.rot / alpha.orth)))
    if (!is.null(sortby)) {
        keys = list()
        for (key in sortby)
            keys[[key]] = with(grid, eval(parse(text=key)))
        ord = do.call("order", keys)
    }
    else
        ord = 1:nrow(grid)
    lim = min(best, nrow(grid))
    ans = grid[ord[1:lim], ]
    row.names(ans) = 1:lim
    ans
}
