KmeansPlus.RNASeq=function(data,nK,model="nbinom",print.steps=FALSE){
    n = data$Count
    s = data$Normalizer
    t = data$Treatment
    logFC = data$logFC
    v = data$NB.Dispersion
    n[n <= 0] = 1e-06
    nG = nrow(n)
    nT = length(unique(t))
    m = sumRow(n, by = t)/sumRow(exp(s), by = t)
    m = m[, t] * exp(s)
    if (model == "poisson") 
        lglk0 = rowSums(-n + m * log(m) - lgamma(n + 1))
    if (model == "nbinom") 
        lglk0 = rowSums(lgamma(n + 1/v) - lgamma(n + 1) - lgamma(1/v) - 
            log(1 + m * v)/v - log(1 + 1/m/v) * n)
    C = matrix(0, nK, nT)
    colnames(C) = colnames(logFC)
    ID = rep(0, nK)
    c0 = rep(0,nT)
    lglk=c()
    for (i in 1:nK) {
        if (print.steps) print(paste("--> centroid",i))
        if (model == "poisson") 
            lglk =cbind(lglk, lglk.ps.c(n = n, s = s, t, C = c0))
        if (model == "nbinom") 
            lglk =cbind(lglk, lglk.nb.c(n = n, s = s, t, v = v, C = c0))
        d = -2 * (lglk - lglk0)
        if (i == 1) 
            d = d[, 1]
        if (i > 1) 
            d = apply(abs(d), 1, min)
        pr = d^2/sum(d^2)
        pr = pr - min(pr) + (max(pr) - min(pr)) * 1e-10
        id = sample(nG, 1, prob = pr)
        C[i, ]= c0 = logFC[id, ]
        ID[i] = id
    }
    class(C) = "logFC"
    return(list(centers = C, id = ID))
}
