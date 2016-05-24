`boot.mlds` <- function(x, nsim, no.warn = TRUE,  ...) {
	if (no.warn){
		old.opt <- options(warn = -1)
		on.exit(old.opt)
		}
	d <- if (x$method == "glm") 
				as.matrix(x$obj$data[, -1]) else
				as.matrix(make.ix.mat(x$data)[, -1])
#				as.mlds.df(x$obj$data) else
#				x$data
	p <- fitted(x)
	rsim <- matrix(rbinom(length(p) * nsim, 1, p), 
					nrow = length(p), ncol = nsim)
	bts.samp <- apply(rsim, 2, function(y, dd) {
#		dd$resp <- x
#		psct <- mlds(dd, ...)$pscale
		
		psct <- glm.fit(dd, y, family = binomial(x$link), ...)$coefficients
		names(psct) <- x$stimulus[-1]
		c(psct, sigma = 1)/psct[length(psct)]
		},
		 dd = d)
	res <- list(boot.samp = bts.samp,
		 bt.mean = apply(bts.samp, 1, mean),
		 bt.sd = apply(bts.samp, 1, sd),
		 N = nsim
		 )
	class(res) <- c("mlds.bt", "list")
	res
	}
	
`boot.mlbs` <- function (x, nsim, no.warn = TRUE, ...) 
{
	if (no.warn){
		old.opt <- options(warn = -1)
		on.exit(old.opt)
	}
    d <- as.matrix(x$obj$data[, -1])
    p <- fitted(x)
    rsim <- matrix(rbinom(length(p) * nsim, 1, p), nrow = length(p), 
        ncol = nsim)
    bts.samp <- apply(rsim, 2, function(y, dd) {
        psct <- glm.fit(dd, y, family = binomial(x$link), ...)$coefficients
        names(psct) <- x$stimulus[-1]
        c(psct, sigma = 1)/psct[length(psct)]
    }, dd = d)
    res <- list(boot.samp = bts.samp, bt.mean = apply(bts.samp, 1, mean), 
        bt.sd = apply(bts.samp, 1, sd), N = nsim)
    class(res) <- c("mlds.bt", "list")
    res
}

`summary.mlds.bt` <- function(object, standard.scale = TRUE, sigma = FALSE, ...) {
	n <- nrow(object$boot.samp)
	if (sigma && standard.scale) 
		print(c(object$bt.mean[n], 
			SD = object$bt.sd[n]))
	if (standard.scale){
		with(object, 
			cbind(mean = c(0, bt.mean[-n]), 
				SD = c(0, bt.sd[-n])))
	} else {
		unnorm <- with(object, 
				apply(boot.samp, 2, function(x) x[-n]/x[n]))
		cbind(means = c(0, rowMeans(unnorm)), 
			SD = c(0, sd(t(unnorm))))
	}
}

