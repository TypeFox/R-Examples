elc <-
function(ax, ag, bx, kt){
    n <- length(kt)   # number of calendar years
    k <- length(ax)   # number of age groups
    g <- length(ag)   # number of extra param. groups
    if (length(bx) > 1) if (length(bx) != k)
        stop('mismatch in age-sp. parameter lengths\n')
    # factor period:
    ft <- gl(n, k)
    fg <- gl(g, n*k)
    array(ax + ag[fg] + bx*kt[ft], dim=c(k, n, g),
          dimnames=list(names(ax), names(kt), names(ag)))
}
