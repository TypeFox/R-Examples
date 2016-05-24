gng.plot.mix <-
function(obj, amplify = 1, resolution=100, new.plot = TRUE,...)
{
	opt <- list(...);
	if(is.null(opt$type))
	{
		type = 'l';
	} else {
		type = opt$type;
	}
	# generate x-points to draw the plot (range of observed data)
	x <- c(seq(obj$range[1],obj$range[2],len=resolution))
	n <- length(x);
	I1 <- (x < (-obj$th1))+0;        
	I2 <- (x > obj$th2)+0;
	f <- matrix(0, n, 1);
	f0 <- matrix (0, n, 1);
	for (k in 2:(obj$K+1)){
		f0 <- f0 + obj$pi[k]*dnorm(x, obj$mu[k-1], obj$sigma[k-1]);		
	}
	e1 <- obj$pi[1]*(I1)*dexp((-x-obj$th1),rate = 1/obj$beta[1]);
	e2 <- obj$pi[obj$K+2]*(I2)*dexp((x-obj$th2),rate = 1/obj$beta[2]);
	f <- (e1 + f0 + e2);
	if(new.plot)
	{	
		if(is.null(opt$xlim)) opt$xlim = range(x);
		if(is.null(opt$ylim)) opt$ylim = range(f * amplify);
		if(is.null(opt$xlab)) opt$xlab = 'x';
		if(is.null(opt$ylab)) opt$ylab = 'Frequency';
		if(is.null(opt$col))  opt$col = 'red';
		if(is.null(opt$main)) opt$main = 'GNG';
		opt$x = 1;
		opt$y = 1;
		opt$type = 'n';
		do.call(plot,opt);
	}
	opt$x <- x;
	opt$y <- f * amplify;
	opt$type =  type;
	do.call(lines,opt);
}

