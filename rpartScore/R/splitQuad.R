splitQuad <-
function(y, wt, x, parms, continuous,sumwt) {
    # Center y
    n <- length(y)
    y <- y- sum(y*wt)/sum(wt)

    if (continuous) {
	
	temp <- cumsum(y*wt)[-n]

	left.wt  <- cumsum(wt)[-n]
	right.wt <- sum(wt) - left.wt
	lmean <- temp/left.wt
	rmean <- -temp/right.wt
	goodness <- 2*(left.wt*lmean^2 + right.wt*rmean^2)/sumwt
	list(goodness= goodness, direction=sign(lmean))
	}
    else {
	# Categorical X variable
	ux <- sort(unique(x))
	wtsum <- tapply(wt, x, sum)
	ysum  <- tapply(y*wt, x, sum)
	means <- ysum/wtsum

	ord <- order(means)
	n <- length(ord)
	temp <- cumsum(ysum[ord])[-n]
	left.wt  <- cumsum(wtsum[ord])[-n]
	right.wt <- sum(wt) - left.wt
	lmean <- temp/left.wt
	rmean <- -temp/right.wt
	list(goodness= 2*(left.wt*lmean^2 + right.wt*rmean^2)/sumwt,
	     direction = ux[ord])
	}
    }

