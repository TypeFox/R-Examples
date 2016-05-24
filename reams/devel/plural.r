#source('subsampleaps.R')
#source("cvic.r")
require(reams)

# Generate table of AIC and CVIC values for best models
plurmod = function(y, X, nfold=length(y), nboot) {
    pmax = ncol(X)+1
	  nnn = length(y)
    cvicobj = cvic(y=y, x=X, pvec=1:pmax, nfold=nfold)
    eicobj = ic.boot(y=y, X=X, nboot = nboot, pvec=1:pmax)
    leapobj = leaps(X, y, method="r2")
    lw = leapobj$which
    term1vec = c()
    
    for (ss in 1:nrow(lw)) 
        term1vec[ss] = nnn * log(crossprod(lsfit(X[ , lw[ss, ]], y)$resid) / nnn)             
    term1vec = c(nnn * log((nnn-1) * var(y) / nnn), term1vec)
    term2vec = eicobj$term2[rep(1:pmax, c(1,table(leapobj$size)))]
    
    edf.mon = cvicobj$edf.mon[rep(1:pmax, c(1,table(leapobj$size)))]

    tabl = matrix(NA, nrow(lw) + 1, ncol(lw) + 4)
    dimnames(tabl) = list(NULL, NULL)
    dimnames(tabl)[[1]] = c(0, dimnames(lw)[[1]]) 
    tabl[1, 1:ncol(lw)] = rep(FALSE, ncol(lw))
    tabl[-1, 1:ncol(lw)] = lw

    dimnames(tabl)[[2]] = c(dimnames(lw)[[2]], 'AIC', 'AICc', 'CVICmon', "EIC")
    ppp = as.numeric(dimnames(tabl)[[1]]) + 1
    tabl[ , ncol(lw)+1] = term1vec + 2 * ppp
    tabl[ , ncol(lw)+2] = term1vec + nnn * (nnn+ppp) / (nnn-ppp-2)
    tabl[ , ncol(lw)+3] = term1vec + nnn*(nnn+edf.mon)/(nnn-edf.mon-2)
    tabl[ , ncol(lw)+4] = term1vec + term2vec
    tabl
}

# Pick out "better-than-nothing" models from AIC/CVIC table
# nmods = maximum # of such models to include
btn = function(modtable, npred, ic='CVICmon', nmods=10, plot=TRUE, labels=1:npred, cex.axis=1, xlab="Predictors", ylab='Criterion value', main=ic) {
	ic.col = modtable[ , ic]
    null.crit = ic.col[1]
    rows.to.incl = ic.col <= min(null.crit, sort(ic.col)[nmods]) 
    if (sum(rows.to.incl)==1) {
    	otpt = matrix(modtable[rows.to.incl, ], 1, ncol(modtable))
    	dimnames(otpt) = list(NULL, colnames(modtable))
    }
    else otpt = modtable[ic.col <= min(null.crit, sort(ic.col)[nmods]), ]
    if (plot) {
     	plot(c(0,npred), range(otpt[ , ic]), type='n', xlab=xlab, ylab=ylab, main=main, xaxt='n')
     	if ("0" %in% rownames(otpt)[1]) abline(h=otpt[1, ic], lty=3, lwd=1.5)
    	axis(1, at=(1:npred)-.5, labels=labels, las=2, cex.axis=cex.axis, tick=FALSE)
    	for (j in 1:nrow(otpt)) for (k in which(otpt[j,1:npred]==1)) lines(c(k-.9,k-.1), rep(otpt[j,ic],2))
    }
    otpt
}    

# Check whether true predictors included more often than false
# Compares CVIC with AIC and AICc
btnsim = function(n, p.all, p.true, R2, nboot, nsim, nmods=7, plot=FALSE) {
	success.ss = success.a1 = success.a2 = rep(0,nsim)
	for (snum in 1:nsim) {
		cat('\nSimulation', snum, '\n')
	    XY = xy(n, p.all, p.true, R2)
        mtable = plurmod(XY$y, XY$X, nboot)
        a1 = btn(mtable, p.all-1, 'AIC', nmods=nmods, plot=plot)
        a2 = btn(mtable, p.all-1, 'AICc', nmods=nmods, plot=plot)
        ss = btn(mtable, p.all-1, nmods=nmods, plot=plot)
        cat('AIC:\n'); print(a1)
        cat('AICc:\n'); print(a2)
        cat('CVICmon:\n'); print(ss)
        sum.a1 = colSums(a1)
        sum.a2 = colSums(a2)
        sum.ss = colSums(ss)
        success.a1[snum] = min(sum.a1[1:(p.true-1)]) > max(sum.a1[p.true:(p.all-1)])
        success.a2[snum] = min(sum.a2[1:(p.true-1)]) > max(sum.a2[p.true:(p.all-1)])
        success.ss[snum] = min(sum.ss[1:(p.true-1)]) > max(sum.ss[p.true:(p.all-1)])
    }
    list(table(success.a1), table(success.a2), table(success.ss))
}
