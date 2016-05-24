plursim <-
function(n, p.all, p.true, R2, nboot, nsim, resample="subsampling", nmods=7) {
	lr = nchar(resample)
	if (substr(resample, 1, lr)==substr("bootstrap", 1, lr)) resam = "boot"
	else if (substr(resample, 1, lr)==substr("subsampling", 1, lr)) resam = "sub"
	else stop("Argument 'resample' must equal 'bootstrap', 'subsampling', or an abbreviation of one of these")

	success.a1 = success.a2 = success.rs = rep(0,nsim)
	for (snum in 1:nsim) {
		cat('\nSimulation', snum, '\n')
	    XY = xy(n, p.all, p.true, R2)
        mtable = scoremods(XY$y, XY$X, nboot, resample=resam)
        a1 = bestmods(mtable, 'AIC', nmods=nmods, plot=FALSE)
        a2 = bestmods(mtable, 'AICc', nmods=nmods, plot=FALSE)
        rs = bestmods(mtable, nmods=nmods, plot=FALSE)
        cat('AIC:\n'); print(a1)
        cat('AICc:\n'); print(a2)
        cat('IC.', resam, ':\n'); print(rs)
        sum.a1 = colSums(a1)
        sum.a2 = colSums(a2)
        sum.rs = colSums(rs)
        
        success.a1[snum] = min(sum.a1[1:(p.true-1)]) > max(sum.a1[p.true:(p.all-1)])
        success.a2[snum] = min(sum.a2[1:(p.true-1)]) > max(sum.a2[p.true:(p.all-1)])
        success.rs[snum] = min(sum.rs[1:(p.true-1)]) > max(sum.rs[p.true:(p.all-1)])
    }
    
    list(AIC=table(success=success.a1), AICc=table(success=success.a2), resampIC=table(success=success.rs))
}

