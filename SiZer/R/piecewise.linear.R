library('boot');

piecewise.linear <- function(x, y, middle=1, 
                             CI=FALSE, bootstrap.samples=1000, 
                             sig.level=.05){
  
  alpha <- piecewise.linear.simple(x,y,middle);
  w <- x-alpha;
  w[w<0] <- 0;
  model <- lm(y~x+w);
  out <- NULL;
  out$change.point <- alpha;
  out$model <- model;
  out$x <- seq(min(x), max(x), length=200);
  w <- out$x-alpha;
  w[w<0] <- 0;
  out$y <- predict(out$model, data.frame(x=out$x, w=w) );
  out$CI <- CI;
  class(out) <- 'PiecewiseLinear';

  # if the user requests confidence intervals for the change point
  # then we'll bootstrap.  Otherwise just return what we have.
  if(CI == FALSE){
    return(out);
  }else{
    data <- data.frame(x=x, y=y);
    # define a function that will return the statistics we wish to bootstrap     
    my.cp <- function(data, index){
      x <- data[index, 1];
      y <- data[index, 2];
      cp <- piecewise.linear.simple(x,y);
      w <- x-cp;
      w[w<0] <- 0;
      model <- lm(y~x+w);
      out <- c(cp, model$coefficients[2],  model$coefficients[3],
                   model$coefficients[2] + model$coefficients[3] );
      return(out);
    }
    boot.result <- boot(data, my.cp, R=bootstrap.samples);
    out$intervals <- apply(boot.result$t, 2, quantile, probs=c(sig.level/2, 1-sig.level/2));
    colnames(out$intervals) <- c('Change.Point', 'Initial.Slope', 
                          'Slope.Change', 'Second.Slope');
    out$CI <- t(out$CI);
    return(out);
  }
}



# This will perform a search for the best changepoint as defined by the 
# maximum likelihood surface.  It uses a search algorithm as opposed 
# to an exhaustive search.  
#
# the middle parameter can be used to narrow the search for a threshold
# by confining the search to the middle data points.  That is, middle=1 and
# all points along the x-axis (between the smallest and largest x values) 
# are possible thresholds.  Middle=.5 confines the search to the middle 50%
# of the range of x values.
piecewise.linear.simple <- function(x, y, middle=1){

  # This is the function that we wish to optimize (over alpha)
  piecewise.linear.likelihood <- function(alpha, x, y){
    N <- length(x);
    w <- (x-alpha);
    w[w<0] <- 0;
    fit <- lm(y ~ x + w);
    Beta <- coefficients(fit);
    Mu <- Beta[1] + Beta[2]*x + Beta[3]*w;    
    SSE <- sum(fit$residuals^2);
    sigma2 <- SSE/N;                    # MLE estimate of sigma^2
    likelihood <- sum( log( dnorm(y, mean=Mu, sd=sqrt(sigma2)) ) );
    return(likelihood);
  }

  r <- range(x);
  offset <- r * (1-middle)/2;
  low <- min(x)  + offset;
  high <- max(x) - offset;
  temp <- optimize(piecewise.linear.likelihood, c(low, high), x=x, y=y, maximum=TRUE);
  return(temp$maximum);
}
  




print.PiecewiseLinear <- function(x, ...){
  print(paste('Threshold alpha:', x$change.point));
  print("");
  print("Model coefficients: Beta[0], Beta[1], Beta[2]");
  print(x$model$coefficients);
  if( x$CI == TRUE){
    print(x$intervals);
  }
}

plot.PiecewiseLinear <- function(x, xlab='X', ylab='Y', ...){
	plot(x$model$model$x, x$model$model$y, xlab=xlab, ylab=ylab, ...);
	lines(x$x, x$y, col='red');
}

predict.PiecewiseLinear <- function( object, x, ...){
	alpha <- object$change.point;
	beta <- object$model$coefficients;
	w <- x - alpha;
	w[w<0] <- 0;
	yhat <- beta[1] + beta[2]*x + beta[3]*w;
	return(yhat);
}

logLik.PiecewiseLinear <- function(object, ...){
	out <- logLik(object$model)
	attr(out,'df') <- 5  # 3 betas, 1 alpha, 1 sigma
	return(out)
}
