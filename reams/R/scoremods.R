scoremods <-
function(y, X, nboot, nfold=length(y), names=NULL) {
    pmax = ncol(X)+1
    nnn = length(y)
    lnames = if (is.null(names)) names(X) else names
    cvicobj = cvic(y=y, X=X, nfold=nfold)
    eicobj = eic(y=y, X=X, nboot=nboot)
    leapobj = leaps(X, y, method="r2", names=lnames)
    lw = leapobj$which
    npred = ncol(lw)
    term1vec = c()
    for (ss in 1:nrow(lw))
        term1vec[ss] = nnn * log(crossprod(lsfit(X[ , lw[ss, ]], y)$resid) / nnn)
    term1vec = c(nnn * log((nnn-1) * var(y) / nnn), term1vec)
    eicpenvec = eicobj$penalty[rep(1:pmax, c(1,table(leapobj$size)))]
    edf.mon = cvicobj$edf.mon[rep(1:pmax, c(1,table(leapobj$size)))]
    tabl = matrix(NA, nrow(lw) + 1, npred + 5)
    dimnames(tabl) = list(NULL, NULL)
    dimnames(tabl)[[1]] = c(0, dimnames(lw)[[1]])
    tabl[1, 1:npred] = rep(FALSE, npred)
    tabl[-1, 1:npred] = lw
    dimnames(tabl)[[2]] = c(dimnames(lw)[[2]], "AIC", "AICc", "BIC", "EIC", "CVICmon")
    ppp = as.numeric(dimnames(tabl)[[1]]) + 1
    tabl[ , npred+1] = term1vec + 2 * ppp
    tabl[ , npred+2] = term1vec + nnn * (nnn+ppp) / (nnn-ppp-2)
    tabl[ , npred+3] = term1vec + log(nnn) * ppp
    tabl[ , npred+4] = term1vec + eicpenvec
    tabl[ , npred+5] = term1vec + nnn*(nnn+edf.mon)/(nnn-edf.mon-2)
    attr(tabl, "npred") = npred
    tabl
}

