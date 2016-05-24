# Helper functions for WQS

#----------------------------------------------------------------------------------------------
# Function to perform preliminary check of the data
#----------------------------------------------------------------------------------------------
check_xyz <- function(x,y,z){
  class_x <- ifelse(class(x)== "matrix" | class(x)== "data.frame", 0, 1)

  n <- length(y)
  n.x <- dim(x)[1]
  dim_x <- ifelse(n != n.x, 1, 0)

  if(is.null(z)){dim_z <- 0}
  else{
    if(class(z) == "matrix") n.z <- dim(z)[1]
    if(class(z) == "vector") n.z <- length(z)[1]
    dim_z <- ifelse(n != n.z, 1, 0)
  }

  return(c(class_x, dim_x, dim_z))
}


#----------------------------------------------------------------------------------------------
# Create ranked data (returns matrix of quantiles)
#----------------------------------------------------------------------------------------------
  quantile.fn <- function(data, n.quantiles){
	q <- matrix(0, dim(data)[1], dim(data)[2])
	I <- dim(data)[2]
	for(i in 1:I){
	q[,i]<- cut(data[,i], breaks = quantile(data[,i], probs = c(0:n.quantiles/n.quantiles)), include.lowest=TRUE)}
	q <- q-1
	#colnames(q) <- paste0("q", seq(1:I))
	return(q)
  }


#----------------------------------------------------------------------------------------------
# objective function for continuous response (returns least squares)
#----------------------------------------------------------------------------------------------
# The objective functions will take the following arguments
# q:	matrix of quantiles (ranked data)
# z:	matrix of covariates
# y:	response vector
# param:	paramter vector (b0,b1,w1,...wc)

# Note: The objective function and the constraint functions must take the same arguments 
#       regardless of whether or not they are used by the particular function

  objfn.cont <- function(param, q, z, y){
	c <- dim(q)[2]
	z.null <- ifelse(is.null(z),1,0) # 1 if NULL
	b0 <- param[1] 	# intercept
	b1 <- param[2] 	# coefficient for WQS term	
	w <- param[3:(2+c)] 	# vector of weights (length c)
	ls <- numeric() # initialize space

	if(z.null == 0){
	  p <- dim(z)[2]	# of covariates
	  theta <- param[(3+c):(2+c+p)] 	# parameters for covariates (length p)
	  mu <- b0 + b1*q%*%w + z%*%theta
	}else{	mu <- b0 + b1*q%*%w}
	ls <- (y-mu)**2
	leastsq <- sum(ls)
	return(leastsq) # minimizes objective fn and we want to minimize least squares
  }


#----------------------------------------------------------------------------------------------
# Linear constraint (allows us to contstrain weights to sum to 1)
#----------------------------------------------------------------------------------------------	
# Note: The objective function and the constraint functions must take the same arguments 
#       regardless of whether or not they are used by the particular function

  lincon <- function(param, q, z, y){
	c <- dim(q)[2]
	weights <- param[3:(2+c)]
	sum <- sum(weights)
	return(sum)
  }


#----------------------------------------------------------------------------------------------
# Inequality constraints (allows us to put constraints on b1 and the weights)
#----------------------------------------------------------------------------------------------	
  ineq <- function(param, q, z, y){
	c <- dim(q)[2]
	b1 <- param[2]
	weights <- param[3:(2+c)]
	return(c(b1, weights))
  }


#--------------------------------------------------------------------------------------------------
# Calculate weights based on relative test statistic
#--------------------------------------------------------------------------------------------------
  teststat.fn <- function(wts, test_stat){
    Sb = abs(test_stat)  # take absolute value so we can calculate relative strength of test statistic
    signal <- Sb/sum(Sb)
    w <- colSums(signal*wts)
    return(w)
  }

#--------------------------------------------------------------------------------------------------
# Function to apply weights to validation data set
#--------------------------------------------------------------------------------------------------
wqs.fit <- function(q, z, y, w){

  WQS <- as.vector(q%*%w)
  z.null <- ifelse(is.null(z),1,0) #1 if NULL
  if(z.null == 0){
  temp <- data.frame(cbind(y, z, WQS))
  } else{temp <- data.frame(cbind(y, WQS))}

  fit <- glm2(y ~ ., data = temp, family = "gaussian"(link = identity))
  out <- list(WQS, fit)
  names(out) <- c("WQS", "fit")
  return(out)
}

#--------------------------------------------------------------------------------------------------
# Function to specify lower and upper bounds
#--------------------------------------------------------------------------------------------------
specify.bounds <- function(b1.pos, c){

  # Specify lower and upper bounds for b1 based on direction of constraint
  if(b1.pos == TRUE){
    b1.LB <- 0 
    b1.UB <- Inf
  } else{
    b1.LB <- -Inf
    b1.UB = 0
  }

  # first term represents bound for b1, next c terms represent bounds for the weights
  ineqLB <- c(b1.LB, rep(0,c))
  ineqUB <- c(b1.UB, rep(1,c))

  out <- list(ineqLB, ineqUB)
  names(out) <- c("ineqLB", "ineqUB")
  return(out)
}

#--------------------------------------------------------------------------------------------------
# Function to specify initial values
#--------------------------------------------------------------------------------------------------
specify.init <- function(z, y, b1.pos, c){

  z.null <- ifelse(is.null(z),1,0) #1 if NULL
  b0.0 <- 0
  w.0 <- rep(1/c, c)
  names(w.0) <- paste0("w", 1:c)
  b1.0 <- ifelse(b1.pos == TRUE, 0.1, -0.1)

  # Initial values for covariates
  if(z.null == 0){
    fit.init <- glm2(y ~ z, family = "gaussian"(link = identity))
    init.z <- coef(fit.init)[-1]
    p <- dim(z)[2]
    names(init.z) <- c(paste0("z",1:p))
    init <- c(b0 = b0.0, b1 = b1.0, w.0, init.z)
  } else{init <- c(b0 = b0.0, b1 = b1.0, w.0)}

  return(init)
}

#--------------------------------------------------------------------------------------------------
# Function to estimate weights across bootstrap samples for WQS Regression
#--------------------------------------------------------------------------------------------------
wqs_b.est <- function(y, q, z, B, pars, fun, eqfun, eqB, ineqfun, ineqLB, ineqUB, LB, UB){ 

	z.null <- ifelse(is.null(z),1,0) #1 if NULL
	c <- dim(q)[2]
	if(z.null == 0){p <- dim(z)[2]} else{p <-0}

      # initialize matrix for parameter estimates (from estimation step)
      result <- matrix(0, nrow = B, ncol = length(pars))
	colnames(result) <- names(pars)
	convergence <- rep(0, B) #0 indicates convergence; 1 or 2 indicates non-convergence

	beta_1 <- rep(0, B)
	pval <- rep(0, B)
	test_stat <- rep(0, B) # test statistic (z-value)

	#---------------------------- BOOTSTRAP ROUTINE -----------------------------------	
	for (b in 1:B) {

	   # draw random sample (of same size as training data set) with replacement
	   samp <- sample(1:length(y), replace = TRUE) 

	   y.b <- as.vector(y[samp])
	   q.b <- q[samp,]   
	   
	   if(z.null == 0){z.b <- as.matrix(z[samp,])}else{z.b <- NULL}
	   rownames(z.b) <- NULL
         result.b <- solnp(pars, fun, eqfun, eqB, ineqfun, ineqLB, ineqUB, LB, UB, q.b, z.b, y.b, 
			control = list(tol = 1e-10,delta = 1e-10, trace = 0))

	   result[b,] <- result.b$pars
	   convergence[b] <- result.b$convergence
         
	   w <- result.b$pars[3:(2+c)]
	   fit <- wqs.fit(q.b, z.b, y.b, w)$fit

         beta_1[b] <- fit$coefficients[row.names = "WQS"]
	   test_stat[b] <- summary(fit)$coefficients['WQS', 't value']  # extract z-value for WQS parameter
	   pval[b] <- summary(fit)$coefficients['WQS', 'Pr(>|t|)'] # extract p-value  
	}

      wts.matrix <- result[,3:(2+c)]; colnames(wts.matrix) <- paste0('w', 1:c)
	out <- list(wts.matrix, convergence, beta_1, test_stat, pval)
	names(out) <- c("wts.matrix", "convergence", "beta_1", "test_stat", "pval")

	return(out)
}
    

