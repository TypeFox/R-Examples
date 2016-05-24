inudge.plot.mix <-
function (obj, amplify = 1, resolution = 100, new.plot = TRUE,...)
{
	opt <- list(...);
	if(is.null(opt$type))
	{
		type = 'l';
	} else {
		type = opt$type;
	}
	# generate x-points to draw the plot (range of observed data)
  x <- c(seq(obj$a, obj$b,len=resolution));
  d <- rep(0,length(x));
  if (abs(obj$a-obj$b) > (1e-15)){
  		d <- (obj$pi[1])*dunif(x,obj$a,obj$b)
  	for (k in 1:obj$K){
  		d <- d + obj$pi[k+1] *dnorm(x, obj$mu[k], obj$sigma[k]);
  	}
  }else {
  	for (k in 1:obj$K){
  		d <- d + obj$pi[k+1] *dnorm(x, obj$mu[k], obj$sigma[k]);
  	}
  }
  	if(new.plot)
	{	
		if(is.null(opt$xlim)) opt$xlim = range(x);
		if(is.null(opt$ylim)) opt$ylim = range(d * amplify);
		if(is.null(opt$xlab)) opt$xlab = 'x';
		if(is.null(opt$ylab)) opt$ylab = 'Frequency';
		if(is.null(opt$col))  opt$col = 'red';
		if(is.null(opt$main)&&obj$name=="iNUDGE") opt$main = 'iNUDGE';
	  if(is.null(opt$main)&&obj$name=="NUDGE") opt$main = 'NUDGE';
		opt$x = 1;
		opt$y = 1;
		opt$type = 'n';
		do.call(plot,opt);
	}
	opt$x <- x;
	opt$y <- d * amplify;
	opt$type =  type;
	do.call(lines,opt);
}

