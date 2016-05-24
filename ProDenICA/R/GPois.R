GPois <-
function(x, df = 6, B = 500, order = 1, widen = 1.2, density.return=FALSE, ...)
{
	### x must be a unit vector with mean zero
	### Note; if x has outliers, this will make the gaussian more robust
	### order and widen are robustness parameters in computing the grid
	### This function also computes derivitives
	x <- drop(scale(x))
	n <- length(x)
	rangex <- range(x)
	if(order == 1)
		rx <- rangex
	else {
		rx <- sort(x)[c(order, n - order + 1)]
	}
	rx <- ylim.scale(rx, diff(rx) * widen)
	xg <- seq(from = rx[1], to = rx[2], length = B)
	gaps <- diff(rx)/(B - 1)
	xcuts <- c(min(rangex[1], rx[1]) - gaps/2, xg[ - B] + gaps/2, max(
		rangex[2], rx[2]) + gaps/2)
	ys <- as.vector(table(cut(x, xcuts)))
	gxg <- dnorm(xg)
	bigdata <- list(xg=xg, gxg=gxg, ys=ys)
	##  assign("bigdata",bigdata,frame=0)
	## assign("df", df, frame = 0)
	pois.fit <- gam(ys ~ s(xg, df) + offset(logb(gxg)), family = poisson,
		data = bigdata, ...)
	## Now to get the derivitives
	Gs <- predict(pois.fit)-logb(gxg)
        if(density.return){
        ## package up the function to return for plotting
        density=list(x=xg,y=exp(Gs+logb(gxg)))
      }
	Gs <- Gs + logb(sum(gxg)/sum(fitted(pois.fit)))
	resp <- Gs + residuals(pois.fit, type = "working")
	weights <- pois.fit$weights
	df <- B - pois.fit$df.residual
	pois.refit <- smooth.spline(xg, resp, weights, df)
	Gs <- predict(pois.refit, x, deriv = 0)$y
	gs <- predict(pois.refit, x, deriv = 1)$y
	gps <- predict(pois.refit, x, deriv = 2)$y
	rl=list(Gs = Gs, gs = gs, gps = gps)
        if(density.return)rl$density=density
        rl
}

