l.actu.data.sp <-
function(l, d, d.sp)
{
    pos <- match(d$loc[l], d.sp$loc)
    z.actu <- d$z[l]
    n.actu <- 1

    return(c(pos, z.actu, n.actu))
}
