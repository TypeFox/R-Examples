powlinear.sel <-
function(imarks, jmarks, dists, dranks, par = list(ki=0.2, kj=0, p=1, r0=0, smark=1)) {
# General competitor selection:  R < ki Si^p + kj Sj^p + r0
    with(as.list(par),
         dists < ki * imarks[[smark]]^p + kj * jmarks[smark]^p + r0
    )
}
