
ExpImprovement<- function(x_new, fmin, fit.gp, muX=NULL, muY=NULL){	

##function EI = ExpImprovement(tst, P)		
	
# fmin = Q_min	

###  x_new - a new point in parameter space

	 require(mlegp)
    # don't allow negative parameters coused by substration 'x_new - muX'!   skip normalisation in R implementation! 
	#x_new= (x_new - muX);

	# check: x_new should be a vector (1, k) !
   # print(" x_new")
    #print((x_new))
	if ((nrow(x_new)) != 1 ) x_new<- t(x_new)
	if ((nrow(x_new)) != 1 ) stop("x_new should be only one point in a k-dim space! dim(x_new) = c(1, k) ") 

    # calculate expected model value E[ y_new | y_obs ] and variance for y_new Var[ y_new | y_obs ]
	pred.list<-predict(fit.gp, x_new, se.fit=TRUE )
    # mu_y_new = mu_y_new + muY 
    mu_y_new <- pred.list[[1]]
    var_y_new <- pred.list[[2]]
    s<- sqrt(abs(var_y_new))
	
    if (abs(fmin - mu_y_new) < 0.01){
		EI <- 0
	}else{
		EI <- (fmin - mu_y_new)*pnorm(fmin, mean = mu_y_new, sd = s) + s*dnorm(fmin, mean = mu_y_new, sd = s)
	}
		EI = -EI; # %% minimization instead of maximization!!!

return(EI)
}


# equal to dnorm!!!!
.normpdf <-function (x,mu,s) {
	#Normal probability density function
	return(	y = 1/(sqrt(2*pi)*s)*exp(-(mu-x)^2/(2*s^2)) )
}


# equal to pnorm !!!!
.normcdf<- function (x,mu,s){
	#Normal cumulative distribution function
	return(y = 0.5  + 0.5 * .erf((x-mu)/(sqrt(2)*s)) )
}

## if you want the so-called 'error function'
.erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1



