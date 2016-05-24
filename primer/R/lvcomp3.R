`lvcomp3` <-
function (t, n, parms) 
{
    with(as.list(parms), {
        dn1dt <- r1 * n[1] * (1 - a11 * n[1] - a12 * n[2] - a13 * 
            n[3])
        dn2dt <- r2 * n[2] * (1 - a22 * n[2] - a21 * n[1] - a23 * 
            n[3])
        dn3dt <- r3 * n[3] * (1 - a33 * n[3] - a31 * n[1] - a32 * 
            n[2])
        list(c(dn1dt, dn2dt, dn3dt))
    })
}
