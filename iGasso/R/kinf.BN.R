kinf.BN = function(G, whole=FALSE)
{
    p = apply(G, 1, mean, na.rm=T)/2                            # missing value in G is allowed
    sel = (p>0.01) & (p<0.99)                                   # only use SNP with MAF > 0.01
    G = G[sel,]
    p = p[sel]

    w = sweep(G, 1, 2*p)    w = sweep(w, 1, sqrt(2*p*(1-p)), FUN="/")
    temp = !is.na(w)    count = crossprod(temp,temp)        # of pairs with non-missing cross product    w[is.na(w)] = 0                     # w'w: missing value does not contribute to the sum
    GG = crossprod(w, w)

    if (!whole) list(total = GG, count = count)      # to be averaged later after all chromosomes are processed
    else GG/count
}
