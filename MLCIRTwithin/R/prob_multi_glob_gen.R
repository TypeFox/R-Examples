prob_multi_glob_gen <- function(X, model, be, ind=(1:dim(X)[3])){
	
# compute global logit probabilities
	nd = dim(X)[3]
    lv = rep(0,nd)
    for(h in 1:nd) lv[h] = sum(!is.na(X[,1,h]))+1
    lm = max(lv)
    ncov = ncol(X)
    nd = dim(X)[3]
    P = matrix(NA, nd, lm)
    for (h in 1:nd) {
	    if (model == "g" || model == "l") {
	    	out = matr_glob(lv[h],model)
	    	Co = out$Co; Ma = out$Ma
	    }
	    G = rbind(matrix(0, 1, lv[h] - 1), diag(lv[h] - 1))
        Xh = matrix(X[1:(lv[h]-1), , h], lv[h] - 1, ncov)
        if (model == "m") {
            P[h, ] = exp(G %*% (Xh %*% be))
            P[h, ] = P[h, ]/sum(P[h, ])
        }
        if (model == "g" | model == "l") 
            P[h,1:lv[h]] = inv_glob(Xh %*% be, model)$p
    }
    Pdis = P; P = Pdis[ind,]
    out = list(Pdis = Pdis, P = P)
    return(out)
    
}