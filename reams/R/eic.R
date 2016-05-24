eic <-
function(y, X, nboot, pvec=1:(ncol(X)+1), say.which=FALSE, reuse=FALSE) {
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
    yboot = matrix(NA, nnn, nboot)
    xboot = array(dim=c(nnn, ncol(X), nboot))
    num = den = matrix(NA, nboot, np)
    for (k in 1:nboot) {
    	samp = sample(nnn, replace=TRUE)
    	yboot[ , k] = y[samp]
        xboot[ , , k] = as.matrix(X[samp, ])        
        lboot = if (reuse) ll else leaps(xboot[ , , k], yboot[ , k], method="r2", nbest=1)
        for (l in 1:np) {
            if (pvec[l]==1) {
            	mod.kl = lm(yboot[ , k] ~ 1)
            	num[k, l] = crossprod(y - mod.kl$coef)
            }
            else {
                if (say.which) cat((1:ncol(X))[lboot$which[lboot$size==pvec[l], ]], '\n')
            	mod.kl = lsfit(xboot[ , lboot$which[lboot$size==pvec[l], ], k], yboot[ , k])
                num[k, l] = crossprod(y - as.matrix(cbind(1, X[ , lboot$which[lboot$size==pvec[l], ]])) %*% mod.kl$coef)
            }
            den[k, l] = crossprod(mod.kl$resid) / nnn
        }
    }
    penalty = colMeans(num / den)
    eic = nlogsig2hat + penalty
    # pstar1 = (penalty*(nnn-2) - nnn^2) / (penalty + nnn)
    names(eic) = # names(pstar1) = 
            names(nlogsig2hat) = names(penalty) = pvec
    list(nlogsig2hat=nlogsig2hat, penalty=penalty, eic=eic, #pstar1=pstar1, 
    best=if (pvec[which.min(eic)]>1) ll$which[ll$size==pvec[which.min(eic)], ] else rep(FALSE, ncol(X)))
}

