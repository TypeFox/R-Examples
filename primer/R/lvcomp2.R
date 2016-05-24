`lvcomp2` <-
function (t, n, parms) 
{
    with(as.list(parms), {
        dn1dt <- r1 * n[1] * (1 - a11 * n[1] - a12 * n[2])
        dn2dt <- r2 * n[2] * (1 - a22 * n[2] - a21 * n[1])
        list(c(dn1dt, dn2dt))
    })
}
