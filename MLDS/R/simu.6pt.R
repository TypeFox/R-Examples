`simu.6pt` <-
function(obj, nsim = 1, nrep, no.warn = TRUE) {
	if (no.warn){
		old.opt <- options(warn = -1)
		on.exit(old.opt)
	}
	dd <- if (obj$method == "glm") 
		as.mlds.df(obj$obj$data) else
		obj$data
	sigma <- obj$sigma
# calculate likelihood by 6-point test
	Six.Pts <- Get6pts(obj, nrep)
	cc <- attr(Six.Pts, "indices")
	ll6pt <- drop(lik6pt(obj, Six.Pts))
# compute p-values for each trial
	pvalues <- fitted(obj)

# 1) rbinom p-values to generate new sets of responses
	rsim <- matrix(rbinom(length(pvalues) * nsim, 1, pvalues), 
			nrow = length(pvalues), ncol = nsim)
# 2) fit glm to each new set (sapply?)
	glc = glm.control(maxit = 50000, epsilon = 1e-14)
	lik.bt <- apply(rsim, 2, function(x, y, d, sp, cc) {
		d$resp <- x
		sp$A$resp <- x[cc[, 1]]
		sp$B$resp <- x[cc[, 2]]
		sp$E$resp <- x[cc[, 3]]
		z <- mlds(d, y$stimulus, lnk = y$link)
		drop(lik6pt(z, sp))
		}, y = obj, d = dd, sp = Six.Pts, cc = cc)

# 3) evaluate original 6-point likelihood w/ respect to 
#     the distribution generated in step 3

boot.lst <- list(boot.samp = lik.bt, lik6pt = ll6pt, 
		p = sum(lik.bt <= ll6pt)/length(lik.bt),
		N = nsim
		) 
boot.lst
}

