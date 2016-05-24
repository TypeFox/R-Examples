ic.min <-
function(y, X, pvec=1:(ncol(X)+1)) {
	require(leaps)
    nnn = length(y)
	np = length(pvec)
    ll = leaps(X, y, method="r2", nbest=1)
    llist = list()
    nlogsig2hat = c()
    for (l in 1:np) {
    	llist[[l]] = if (pvec[l]==1) lsfit(rep(1,nnn), y, intercept=FALSE) else lsfit(X[ , ll$which[ll$size==pvec[l], ]], y)
        nlogsig2hat[l] = nnn * log(crossprod(llist[[l]]$resid) / nnn)
    }
    aic = nlogsig2hat + 2 * pvec
    aicc = nlogsig2hat + nnn * (nnn + pvec) / (nnn - pvec - 2)
    bic = nlogsig2hat + pvec*log(nnn)
    list(nlogsig2hat=nlogsig2hat, aic=aic, aicc=aicc, bic=bic,
         best.aic=if (pvec[which.min(aic)]>1) ll$which[ll$size==pvec[which.min(aic)], ] else rep(FALSE, ncol(X)),
         best.aicc=if (pvec[which.min(aicc)]>1) ll$which[ll$size==pvec[which.min(aicc)], ] else rep(FALSE, ncol(X)),
         best.bic=if (pvec[which.min(bic)]>1) ll$which[ll$size==pvec[which.min(bic)], ] else rep(FALSE, ncol(X)))
}

