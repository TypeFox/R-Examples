cvic <-
function(y, X, nfold=length(y), pvec=1:(ncol(X)+1))   {
    nnn = length(y)
	if (nnn %% nfold != 0) stop("Sample size must be divisible by 'nfold'")
    require(leaps)
    all.folds = split(sample(nnn),  rep(1:nfold, length=nnn))
    np = length(pvec)
    n.train = nnn - nnn/nfold
    llist = list()
    nlogsig2hat = c()
    num = den = matrix(NA, nrow=nfold, ncol=np)
    ll = leaps(X, y, method="r2", nbest=1)
    for (l in 1:np) {
	      llist[[l]] = if (pvec[l]==1) lsfit(rep(1,nnn), y, intercept=FALSE) else lsfit(X[ , ll$which[ll$size==pvec[l], ]], y)
        nlogsig2hat[l] = nnn * log(crossprod(llist[[l]]$resid) / nnn)
    }
    for (k in 1:nfold)    {
        idx = all.folds[[k]]
        ltrain = leaps(X[-idx,], y[-idx], method="r2", nbest=1)
        for (l in 1:np)   {
            if (pvec[l]==1) {
            	mod.kl = lm(y[-idx] ~ 1)
            	num[k, l] = crossprod(y[idx] - mod.kl$coef)
            }
            else {
                mod.kl = lsfit(X[-idx, ltrain$which[ltrain$size==pvec[l], ]], y[-idx])
                xval = cbind(1, matrix(X[idx, ltrain$which[ltrain$size==pvec[l], ]], nrow=length(idx)))
                if (length(idx)==1) xval = as.numeric(xval)
                num[k, l] = crossprod(y[idx] - xval %*% mod.kl$coef)
            }
            den[k, l] = crossprod(mod.kl$resid) / n.train
        }
    }
    cv.pen = colSums(num / den) 
    # pen2df = function(pen) min(Re(polyroot(c(pen*(n.train-1)*(n.train-2)-(n.train-1)*n.train*nnn, pen*(3-2*n.train), pen))))
    pen2df = function(pen) n.train - sqrt(nnn*(n.train+1)*(n.train-2) / pen) - 2
    edf = sapply(cv.pen, pen2df)
    
    # Constrained monotone smoothing of edf
    require(mgcv)
    nk = min(30,np)
    dat<-data.frame(pvec=pvec, edf=edf)
    f.ug<-gam(edf~s(pvec,k=nk,bs="cr"),method="REML")
    sm<-smoothCon(s(pvec,k=nk,bs="cr"),dat,knots=NULL)[[1]]
    F<-mono.con(sm$xp)   # monotonicity constraints
    G<-list(X=sm$X,C=matrix(c(1,rep(0,2*(nk-1)),1), 2, nk),sp=f.ug$sp,p=sm$xp,y=edf,w=edf*0+1)
    G$Ain<-F$A;G$bin<-F$b;G$S<-sm$S;G$off<-0
    p<-pcls(G);  # fit spline (using s.p. from unconstrained fit)
    edf.mon <- as.vector(Predict.matrix(sm,data.frame(pvec=pvec))%*%p)
    
    cvic = nlogsig2hat + nnn*(nnn+edf)/(nnn-edf-2)
    cvic.mon = nlogsig2hat + nnn*(nnn+edf.mon)/(nnn-edf.mon-2)
    names(cvic) = names(cvic.mon) = names(edf) = names(edf.mon) = names(nlogsig2hat) = names(cv.pen) = pvec
    list(nlogsig2hat=nlogsig2hat, cv.pen=cv.pen, edf=edf, edf.mon=edf.mon, cvic=cvic, cvic.mon=cvic.mon, 
    best=if (pvec[which.min(cvic)]>1) ll$which[ll$size==pvec[which.min(cvic)], ] else rep(FALSE, ncol(X)),
    best.mon=if (pvec[which.min(cvic.mon)]>1) ll$which[ll$size==pvec[which.min(cvic.mon)], ] else rep(FALSE, ncol(X)))
}

