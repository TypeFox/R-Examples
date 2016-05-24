bestmods <-
function(tbl, ic="AIC", nmods=10, plot=TRUE, 
labels=dimnames(tbl)[[2]][1:attr(tbl,"npred")], cex.axis=1, 
las=2-all(labels==1:attr(tbl,"npred")), 
xlab=if (las==1) "Predictors" else "", ylab="Criterion value", main=ic, ...) {
	nmods = min(nmods, nrow(tbl))
	npred = attr(tbl, "npred")
	ic.col = tbl[ , ic]
    null.crit = ic.col[1]
    rows.to.incl = ic.col <= min(null.crit, sort(ic.col)[nmods])
    if (sum(rows.to.incl)==1) {
    	otpt = matrix(tbl[rows.to.incl, ], 1, ncol(tbl))
    	dimnames(otpt) = list(NULL, colnames(tbl))
    }
    else otpt = tbl[ic.col <= min(null.crit, sort(ic.col)[nmods]), ]
    if (plot) {
     	plot(c(0,npred), range(otpt[ , ic]), type="n", xlab=xlab, ylab=ylab, main=main, xaxt="n", ...)
     	if ("0" %in% rownames(otpt)[1]) abline(h=otpt[1, ic], lty=3, lwd=1.5)
    	axis(1, at=(1:npred)-.5, labels=labels, las=las, cex.axis=cex.axis, tick=FALSE)
    	for (j in 1:nrow(otpt)) for (k in which(otpt[j,1:npred]==1)) lines(c(k-.9,k-.1), rep(otpt[j,ic],2))
    }
    otpt
}

