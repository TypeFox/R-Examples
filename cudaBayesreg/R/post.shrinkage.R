#
# visualizing shrinkage as a function of mean beta{vreg} estimates
#
post.shrinkage.mean <-
function (out, X, vreg, plot=T) 
{
    pmeans <- pmeans.hcoef(out$betadraw)
    yrec <- X %*% (t(pmeans))
    yrecmean <- apply(yrec, 2, mean)
    rg <- range(yrecmean)
    beta <- pmeans[, vreg]
    xlim <- range(beta)
		if(plot)
	    plot(beta, yrecmean, ty = "p", xlim = xlim,
				xlab = paste("beta", vreg), ylab = "post mean", cex = 0.6)
    invisible(list(yrecmean=yrecmean,beta=beta))
}

post.shrinkage.minmax <-
function (out, X, vreg, plot=T) 
{
    pmeans <- pmeans.hcoef(out$betadraw)
    yrec <- X %*% (t(pmeans))
    beta <- pmeans[, vreg]
    xlim <- range(beta)
    yrecmax <- apply(yrec, 2, max)
    yrecmin <- apply(yrec, 2, min)
    ylim <- range(yrecmin, yrecmax)
		if(plot) {
    	plot(beta, yrecmax, ty = "p", ylim = ylim, xlim = xlim, 
        xlab = paste("beta", vreg), ylab = "post max. and min.", 
        pch = 1, cex = 0.6, col = "blue")
    	points(beta, yrecmin, pch = 6, cex = 0.6, col = "red")
		}
    invisible(list(yrecmin=yrecmin, yrecmax=yrecmax, beta=beta))
}
