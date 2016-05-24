#library('splines');

#######################################################################
#  locally.weighted.polynomial
#######################################################################
# Computes the locally weighted polynomial smoothing of the given data
#######################################################################
# Inputs:
# x, y        - bivariate data vectors to be smoothed
# x.grid      - either NA, an integer, or a vector.  If NA, then a default number (41)
#               of grid points is used, spread evenly across the range of the x-vector.
#               If an integer, it designates how many grid points are used, again spread
#               evenly across the range of the x-vector.  If x.grid is a vector, we just
#               use what is given
# h           - The smoothing bandwidth to use, in terms of standard deviations
# degree      - the degree of polynomial that is fit locally, the default is to fit
#               a linear function, ie degree=1
# kernel.type - The shape of the kernal function.  Acceptable values are, 'Normal',
#               'Epanechnikov', 'biweight', and 'triweight'
#######################################################################
# Outputs:
# data   - a structure of the data used to generate the smoothing curve
# h      - the bandwidth used to generate the smoothing curve
# x.grid - the grid of x-values that we have estimated function value
#          and derivative(s) for
# degrees.freedom - the effective sample size at each grid point
# Beta   - a matrix of estimated beta values.  The number of rows is
#          degrees+1, while the number of columns is the same as the length
#          of x.grid.
#          Notice that hat{f}(x_i) = Beta[1,i]
#                      hat{f'}(x_i) = Beta[2,i]*1!
#                      hat{f''}(x_i) = Beta[3,i]*2!  and so on...
# Beta.var - matrix of estimated variances for Beta.  Same structure dimensions.
#############################################################################
locally.weighted.polynomial <- function(x, y, h=NA, x.grid=NA, degree=1, kernel.type='Normal'){

  # convert whatever we were sent in to a valid grid.  If sent a valid grid, we'll keep it
  x.grid <- x.grid.create(x.grid, x,y);
  
  if(is.na(h)){
    h <- bw.SJ(x);
  }
  
  # set up the output data structure
  out <- NULL;
  out$data <- NULL;
  out$data$x <- x;              # storing the data used to calculate this object
  out$data$y <- y;
  out$data$n <- length(x);
  out$h <- h;                   # bandwidth parameter
  out$x.grid <- x.grid;         # where the smoothed function value is to be calculated


  # the effective degrees of freedom for each of the local regressions 
  out$effective.sample.sizes <- rep(NA, length(x.grid));   

  out$Beta <- matrix(nrow=degree+1, ncol=length(x.grid));  # 1st row is function value, beta[0]
                                                           # 2nd row is slope, beta[1]
                                                           # 3rd row is beta[2]  : 
                                                           #    2nd derivative = beta[2] * 2!
											                               #    3rd derivative = beta[3] * 3!
												
	unscaled.var  <- out$Beta;	# Calculating the Variance for the Beta vector will be a 2 step process
	out$Beta.var  <- out$Beta;	# We'll store part while finding f.hat(Y|X), and then use the 
												# residuals to get Var(Beta|X)
	var.y.given.x <- rep(0, length(x.grid) ); 

  # index denoting what where we are in the array x.grid
  count <- 1;

  for(x.star in x.grid){
     # figure out the X matrix
		X <- NULL;
		for( i in 0:degree ){
			X <- cbind(X, (x-x.star)^i );
		}
     
     # find the weight vector
     w <- kernel.h(x-x.star, h, type=kernel.type);
     w2 <- w^2;
     XtW  <- t( apply(X,2,function(x){x*w}) );           # XtW <- t(X) %*% diag(w)   computed faster.  :) 
     XtW2 <- t( apply(X,2,function(x){x*w2}) );          # XtW2 <- t(X) %*% diag(w^2) 
     XtWXinv <- try( solve( XtW %*% X ), silent=TRUE );
     if( class(XtWXinv) == 'try-error' ){            # inverse failed!  Could be any of a 
     		out$Beta[, count] <- rep(NA, degree+1);      # bunch of legitimate reasons such 
        unscaled.var[, count] <- rep(NA, degree+1);  # as a really small bandwidth.  So just 
                                                     # record NAs and don't throw any warnings or errors.   
     }else{
     		beta <- XtWXinv %*% XtW %*% y;
        out$Beta[,count] <- beta;
        unscaled.var[,count] <- diag(XtWXinv %*% XtW2%*%X %*% t(XtWXinv));         
                                                 # Var(Beta_j) = ((XtWX)^{-1}(XtW2X)(XtWX)^-1)[j,j] * 
                                                 #                Var(Y|X=x.star) 
                                                 # but we don't know
     }                                           # Var(Y|X) yet so just store the diagonal terms.
     count <- count+1;
  } 

	# At this point we have calculated the Beta coefficients at each grid point.  Now we need the 
	# variance estimates.
  
	# first get the residuals by interpolating along the grid points
	index <- which(!is.na(out$Beta[1,]))
	if(length(index) > 2 ){
		x.grid.small <- x.grid[index]
		interp <- interpSpline(x.grid.small, out$Beta[1,index])
		out$residuals <- y - predict(interp, x)$y
    
		# Var(Y|X=x) = ( sum_{j=1}^n  e_j^2 * K_h(x-X_j) )    /   
		#        ( sum_{j=1}^n K_h(x-X_j) )  from Chaudhuri and Marron 1999
		for(i in 1:length(x.grid.small) ){
			divisor <- 0
  			for(j in 1:length(x)){
				kernel <- kernel.h(x.grid.small[i]-x[j], h, type=kernel.type)
				var.y.given.x[i] <- var.y.given.x[i] + (out$residuals[j])^2 * kernel
				divisor <- divisor + kernel
 			}
			var.y.given.x[i] <- var.y.given.x[i] / divisor;
		}
	}
	index <- which(var.y.given.x == 0)
	var.y.given.x[index] <- NA
	
	
	# calculate the effective sample size
	for( i in 1:length(x.grid) ){
		sum <- 0;
		for( j in 1:length(x) ){
			sum <- sum + kernel.h(x.grid[i] - x[j], h, type=kernel.type);
		}
		out$effective.sample.sizes[i] <- sum / kernel.h(0, h, type=kernel.type);
	}	 
  
	# now we know everything for calculating the variance of the Beta values
	for( i in 1:length(x.grid) ){
		out$Beta.var[, i] <- unscaled.var[,i] * var.y.given.x[i];
	}
  
  class(out) <- 'LocallyWeightedPolynomial';       # set the out object to have the right label
  return(out);
} 



plot.LocallyWeightedPolynomial <- function(x, derv=0, CI.method=2, alpha=.05, use.ess=TRUE, draw.points=TRUE, ...){
   index <- derv+1;
   intervals <- calc.CI.LocallyWeightedPolynomial(x, derv=derv, CI.method, alpha=alpha, use.ess);
   
   if( any( is.na(intervals$lower.CI) ) ){
   	  stop('Too few data points to perform a valid local regression.  Set use.ess=F to ignore the issue or select a larger bandwidth.');
   }
   # Figure out decent ranges for the plot; User can override these using passing in 
   # vectors for xlim and and ylim 
   
   y.M <- max( intervals$upper.CI );
   y.m <- min( intervals$lower.CI );
   if(derv==0 & draw.points==TRUE){
   	  y.M <- max( c(y.M, x$data$y) );  # If we have to draw the data points, we'll want
   	  y.m <- min( c(y.m, x$data$y) );  # the y.Max and y.min set up to show all the points
   }
   temp <- range(x$x.grid);
   x.m <- temp[1];
   x.M <- temp[2];

   # set up the plot device to be the right size
   plot( c(x.m, x.M), c(y.m, y.M), type='n', ...);

   # draw the confidence interval lines
   lines(x$x.grid, intervals$upper.CI );
   lines(x$x.grid, intervals$lower.CI );
   
   # Shade in between the CI bands
   polygon(c(x$x.grid, rev(x$x.grid)), c(intervals$upper.CI, rev(intervals$lower.CI)),
           col='light grey');

   # draw the middle line (do it after shading otherwise shading is on top)
   lines(x$x.grid, intervals$estimated, type='l', lwd=2);

   # If we are plotting a derivative, then put a horizontal line at zero
   if(derv >= 1){
   	 lines( range(x$x.grid), c(0,0) );
   }

   # if we should display the data points
   if( derv==0 & draw.points == TRUE ){
      points(x$data$x, x$data$y, ...);
   }

}




calc.CI.LocallyWeightedPolynomial <- function(model, derv=0, CI.method=2, alpha=.05, use.ess=TRUE){
   out <- NULL;
   out$lower.CI <- rep(NA, length(model$x));
   out$estimated <- out$lower.CI;
   out$upper.CI <- out$lower.CI;
   
   index <- derv+1;   # which row of the beta matrix we need to deal with

   # Pointwise intervals
   if( CI.method == 1 ){
	   t.star <- qt(1-alpha/2, df=model$degrees.freedom);   # this is a vector!
       out$estimated <- model$Beta[index,];
       out$lower.CI <- out$estimated - t.star*sqrt( model$Beta.var[index,] );
       out$upper.CI <- out$estimated + t.star*sqrt( model$Beta.var[index,] );           
   }
   # Using Hannig and Marron 2006 row wise confidence intervals (eq 10 and 11) 
   if( CI.method == 2 ){
   	   g <- length(model$x.grid);
   	   h <- model$h;
   	   delta <- model$x.grid[2] - model$x.grid[1];
   	   theta <- 2* pnorm( sqrt((1+2*derv)*log(g)) * delta/(2*h) ) -1; 
   	   temp.star <- qnorm( (1-alpha/2)^(1/(theta*g)) );
       out$estimated <- model$Beta[index,] * factorial(derv);
       out$lower.CI <- out$estimated - temp.star*sqrt(model$Beta.var[index,]) * factorial(derv);
       out$upper.CI <- out$estimated + temp.star*sqrt(model$Beta.var[index,]) * factorial(derv);
   }		

   # Put in NAs when the effective sample size is too small for good results
   if( use.ess == TRUE){
      index <- which( model$effective.sample.sizes < 5 );
      out$lower.CI[index] <- NA;
      out$upper.CI[index] <- NA;
   }
   return(out);	
}	


kernel.h <- function(x, h, type='Normal'){
	if(type == 'Normal'){
		out <- dnorm(x/h) / h;
	}else if(type == 'Epanechnikov'){
        out <- (1/beta(.5,2)) * ( positive.part( (1-(x/h)^2) ) );
	}else if( type == 'biweight' ){
		out <- (1/beta(.5,3)) * ( positive.part( (1-(x/h)^2) ) )^2;	}else if( type == 'triweight' ){
		out <- (1/beta(.5,4)) * ( positive.part( (1-(x/h)^2) ) )^3;	}else{   #uniform
		out <- rep(1, length(x));
		out[ abs(x) > h ] <- 0;
	}
    return(out);
} 			



x.grid.create <- function(x.grid, x, y, grid.length=41){
  if(all(is.na(x.grid))){
     # if we are passed an NA
  	 out <- seq(min(x), max(x), length=grid.length);
  }else if( length(x.grid) == 1 ){
     # if we are passed a single integer
     out <- seq(min(x), max(x), length=x.grid);    
  }else{
  	 # otherwise just return what we are sent
  	 out <- x.grid;
  } 
  return(out);
}



	
# for each point along the x-grid, figure out what state we are in
find.states <- function(intervals){
   n <- length(intervals$lower.CI);
   out <- rep(NA, n);
   for( i in 1:n ){
   	  out[i] <- find.state( intervals$lower.CI[i], intervals$upper.CI[i] );
   }
   return(out);
}   
   

# returns the indices of where the state changes occur.  Also returns
# the states that were passed through
find.state.changes <- function(intervals){
    out <- NULL;
    out$indices <- NULL;
    out$state <- NULL;

    lower.CI <- intervals$lower.CI;
    upper.CI <- intervals$upper.CI;

    count <- 1;
    out$state[count] <- find.state(lower.CI[1], upper.CI[1]);
	out$indices[count] <- 1;
	continue <- TRUE;
	
	while(continue == TRUE){
	   if(out$state[count] == 1){
          compare <- lower.CI < 0;
       }else if(out$state[count] == -1){
       	  compare <- upper.CI > 0;
       } else{
       	  compare <- upper.CI < 0 | lower.CI > 0
       }
	   index <- match( TRUE, compare );
	   if( is.na(index) ){
	   	  continue <- FALSE;
	   }else{
          count <- count + 1;
          lower.CI <- lower.CI[ -1*1:(index-1) ];
          upper.CI <- upper.CI[ -1*1:(index-1) ];
          out$state[count] <- find.state(lower.CI[1], upper.CI[1]);           	      out$indices[count] <- out$indices[count-1] + index - 1;
       }
    }
    return(out);   	   	  		
}          	   	   	  		

# given two scalars, find out if the are both positive, one of each, or both negative
# The states are coded as positive=1, negative=-1, possibly zero=0
#  the state 2 corresponds to NotEnoughData, i.e. the effective sample size is < 5
find.state <- function(lower.CI, upper.CI){
    if( is.na(lower.CI) ){
    	state <- 2;
    }else if( lower.CI > 0 ){
		state <- 1;        
	}else if( upper.CI < 0 ){
		state <- -1;
	}else{
		state <- 0;
	}
	
    return(state);
}     	   	   	  		


positive.part <- function(x){
	out <- x;
	out[ out<0 ] <- 0;
	return(out);
}		



