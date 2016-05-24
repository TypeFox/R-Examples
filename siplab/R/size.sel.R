size.sel <-
function(imarks, jmarks, dists, dranks, par = list(k=0.2, smark=1)) {
# Competing radius is k times the subject plant size
    dists < par['k'] * imarks[[par['smark']]]
}
