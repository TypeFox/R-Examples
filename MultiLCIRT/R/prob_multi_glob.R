prob_multi_glob <- function(X, model, be, ind=(1:dim(X)[3])){
	
# compute global logit probabilities	
    l = nrow(X) + 1
    ncov = ncol(X)
    nd = dim(X)[3]
    if (model == "g" || model == "l") {
        out = matr_glob(l, model)
        Co = out$Co
        Ma = out$Ma
    }
    G = rbind(matrix(0, 1, l - 1), diag(l - 1))
    P = matrix(0, nd, l)
    for (h in 1:nd) {
        Xh = matrix(X[, , h], l - 1, ncov)
        if (model == "m") {
            P[h, ] = exp(G %*% (Xh %*% be))
            P[h, ] = P[h, ]/sum(P[h, ])
        }
        if (model == "g" | model == "l") 
            P[h, ] = inv_glob(Xh %*% be, model)$p
    }
    Pdis = P; P = Pdis[ind,]
    out = list(Pdis = Pdis, P = P)
    return(out)
    
}