powers.ker <-
function(imarks, jmarks, dists, dranks, par = list(pi=1, pj=1, pr=1, smark=1)) {
# A generalized competition kernel, (Sj^pj / Si^pi) / dist^pr
    with(as.list(par),
         (jmarks[smark]^pj / imarks[[smark]]^pi) / dists^pr
    )
}
