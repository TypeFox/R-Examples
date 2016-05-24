`plotRate` <-
function(bt, pars)
{
 	time <- seq(0, max(bt), length.out=100);
 	lam <- sapply(time, lambdaFx, lam0=pars$lam0, k = pars$k);
 	mu <- sapply(time, muFx, mu0 = pars$mu0, z = pars$z);
 	plot.new();
 	ymin <- min(lam, mu);
 	ymax <- max(lam, mu);		
	plot.window(xlim=c(0, max(time)), ylim=c(ymin, ymax));
	lines(x=time, y=lam, col='red');
	lines(x=time, y=mu, col='black');
	axis(1);
	axis(2);
}

