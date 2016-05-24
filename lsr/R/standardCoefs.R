# file:    standardCoefs.R 
# author:  Dan Navarro
# contact: daniel.navarro@adelaide.edu.au
# changed: 13 November 2013

# standardCoefs() computes the standardised regression coefficients (i.e. beta 
# weights) for a regression model. It takes an lm object as input, and outputs a
# matrix showing both the raw "b" weights and the standardised "beta" weights. 
# Note that, if you're serious about computing the relative importance of the
# predictors, you'd be better off looking at the relaimpo package, which provides
# better tools for relative importance regression. The standardCoefs() function 
# is just a handy crutch for beginners: beta weights aren't the best solution to
# the problem.
standardCoefs <- function( x ) {
  
  if( !is(x,"lm") ) {stop( '"x" must be a linear model object')}
  
  # read off the useful info
  term.names <- names(x$coefficients)[-1] # all names except "intercept"
  b <- x$coefficients[-1] # grab coefficients except "intercept"
  
  # construct the design matrix
  predictors <- model.matrix(x$terms, data = x$model)
  predictors <- predictors[ ,-1, drop=FALSE] # (we don't want the intercept)
  predictors <- as.data.frame(predictors)  # hack!!
  
  # standard deviations
  sy <- sd(x$model[[1]])       # sd for the outcome  
  sx <- sapply(predictors,sd)  # sd for the predictors
  
  # now compute beta
  beta <- b * sx / sy
  
  # convert to matrix
  coefficients <- cbind(b,beta)
  colnames(coefficients) <- c("b","beta")
  rownames(coefficients) <- term.names
  
  return(coefficients)
  
  
}
